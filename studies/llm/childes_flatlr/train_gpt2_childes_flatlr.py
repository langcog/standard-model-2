"""Train GPT-2 small on CHILDES with a configurable LR schedule, logging val_bpb +
CDI-word surprisal + BLiMP/Zorro grammar accuracy along a shared log-spaced step axis.

Patched fork of feng_eval/train_gpt2_childes.py for the flat-vs-decay LR experiment:
  - --lr_scheduler_type {linear,constant,cosine,...}  (orig hardcoded "linear")
  - grammar trajectory callback (Zorro + BLiMP) at log-spaced steps
  - early stopping OFF by default (we want the full trajectory, not the val-optimal stop)
Deploy in feng_eval/ so surprisal_callback / grammar_callback / eval_grammar import cleanly.
"""
import argparse, os, sys
import torch
from datasets import load_dataset
from transformers import (
    AutoTokenizer, GPT2Config, GPT2LMHeadModel, DataCollatorForLanguageModeling,
    Trainer, TrainingArguments, EarlyStoppingCallback, set_seed,
)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from surprisal_callback import WordSurprisalCallback
from grammar_callback import GrammarTrajectoryCallback


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_file", required=True)
    ap.add_argument("--val_file", required=True)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--config_file", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--cdi_contexts_jsonl", required=True)
    ap.add_argument("--surprisal_csv", required=True)
    ap.add_argument("--num_train_epochs", type=float, default=40.0)
    ap.add_argument("--per_device_batch_size", type=int, default=8)
    ap.add_argument("--learning_rate", type=float, default=1e-4)
    ap.add_argument("--lr_scheduler_type", default="linear",
                    help="linear (warmup+decay-to-0) | constant (FLAT) | cosine | ...")
    ap.add_argument("--warmup_steps", type=int, default=0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n_log_points", type=int, default=60)
    ap.add_argument("--save_total_limit", type=int, default=1)
    ap.add_argument("--seq_len", type=int, default=1024)
    ap.add_argument("--logging_steps", type=int, default=50)
    ap.add_argument("--eval_callback_batch_size", type=int, default=128)
    ap.add_argument("--eval_max_ctx", type=int, default=128)
    ap.add_argument("--eval_max_per_word", type=int, default=50)
    ap.add_argument("--no_save", action="store_true")
    ap.add_argument("--early_stopping_patience", type=int, default=0)
    ap.add_argument("--max_eval_blocks", type=int, default=500)
    # ---- grammar trajectory ----
    ap.add_argument("--zorro_dir", required=True)
    ap.add_argument("--grammar_csv", required=True)
    ap.add_argument("--n_grammar_points", type=int, default=12)
    ap.add_argument("--grammar_max_pairs", type=int, default=100)
    ap.add_argument("--no_blimp", action="store_true")
    ap.add_argument("--grammar_val_blocks", type=int, default=200,
                    help="val LM blocks subsampled for the aligned val_bpb readout")
    args = ap.parse_args()

    set_seed(args.seed)
    print(f"[train] tokenizer {args.tokenizer_dir}", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(args.tokenizer_dir)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    config = GPT2Config.from_pretrained(args.config_file)
    config.vocab_size = tokenizer.vocab_size
    config.bos_token_id = tokenizer.bos_token_id if tokenizer.bos_token_id is not None else tokenizer.eos_token_id
    config.eos_token_id = tokenizer.eos_token_id
    model = GPT2LMHeadModel(config)
    print(f"[train] params: {sum(p.numel() for p in model.parameters()):,} "
          f"(vocab {config.vocab_size}, n_layer {config.n_layer}, n_embd {config.n_embd})", flush=True)

    raw = load_dataset("text", data_files={"train": args.train_file, "validation": args.val_file},
                       keep_linebreaks=True)

    def tok_fn(ex):
        return tokenizer(ex["text"])
    tokenized = raw.map(tok_fn, batched=True, remove_columns=["text"], desc="Tokenizing", num_proc=4)

    block_size = args.seq_len

    def group_texts(ex):
        concat = {k: sum(ex[k], []) for k in ex.keys()}
        total = (len(concat[list(ex.keys())[0]]) // block_size) * block_size
        res = {k: [t[i:i + block_size] for i in range(0, total, block_size)] for k, t in concat.items()}
        res["labels"] = res["input_ids"].copy()
        return res
    lm = tokenized.map(group_texts, batched=True, num_proc=4, desc=f"Grouping {block_size}")
    print(f"[train] blocks train {len(lm['train'])} val {len(lm['validation'])}", flush=True)

    eval_dataset = lm["validation"]
    if args.max_eval_blocks > 0 and len(eval_dataset) > args.max_eval_blocks:
        eval_dataset = eval_dataset.select(range(args.max_eval_blocks))

    # val subsample tensor for the aligned val_bpb readout in the grammar callback
    nvb = min(args.grammar_val_blocks, len(lm["validation"]))
    val_ids = torch.tensor([lm["validation"][i]["input_ids"] for i in range(nvb)], dtype=torch.long)

    early_stop = args.early_stopping_patience > 0
    save_strategy = "epoch" if early_stop else ("no" if args.no_save else "epoch")
    targs = TrainingArguments(
        output_dir=args.output_dir, overwrite_output_dir=True, do_train=True, do_eval=True,
        num_train_epochs=args.num_train_epochs,
        per_device_train_batch_size=args.per_device_batch_size,
        per_device_eval_batch_size=args.per_device_batch_size,
        learning_rate=args.learning_rate,
        lr_scheduler_type=args.lr_scheduler_type,      # <-- the experimental knob
        warmup_steps=args.warmup_steps,
        weight_decay=0.0, adam_beta1=0.9, adam_beta2=0.999, adam_epsilon=1e-8,
        seed=args.seed, data_seed=args.seed, logging_steps=args.logging_steps,
        eval_strategy="epoch", save_strategy=save_strategy, save_total_limit=args.save_total_limit,
        load_best_model_at_end=early_stop,
        metric_for_best_model="eval_loss" if early_stop else None,
        greater_is_better=False if early_stop else None,
        report_to=[], bf16=torch.cuda.is_bf16_supported(), fp16=not torch.cuda.is_bf16_supported(),
        dataloader_num_workers=2, ddp_find_unused_parameters=False,
    )

    class NoShuffleTrainer(Trainer):
        def _get_train_sampler(self):
            return torch.utils.data.SequentialSampler(self.train_dataset)

    cdi_cb = WordSurprisalCallback(
        contexts_jsonl=args.cdi_contexts_jsonl, out_csv=args.surprisal_csv,
        n_points=args.n_log_points, eval_batch_size=args.eval_callback_batch_size,
        max_ctx=args.eval_max_ctx, max_per_word=args.eval_max_per_word)
    gram_cb = GrammarTrajectoryCallback(
        tokenizer=tokenizer, zorro_dir=args.zorro_dir, out_csv=args.grammar_csv,
        n_points=args.n_grammar_points, max_pairs=args.grammar_max_pairs,
        do_blimp=not args.no_blimp, eval_batch_size=args.eval_callback_batch_size,
        val_input_ids=val_ids)

    callbacks = [cdi_cb, gram_cb]
    if early_stop:
        callbacks.append(EarlyStoppingCallback(early_stopping_patience=args.early_stopping_patience))

    trainer = NoShuffleTrainer(
        model=model, args=targs, train_dataset=lm["train"], eval_dataset=eval_dataset,
        tokenizer=tokenizer, data_collator=DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False),
        callbacks=callbacks)
    trainer.train()
    if not args.no_save and not early_stop:
        trainer.save_model(os.path.join(args.output_dir, "final"))
    print(f"[train] DONE. cdi={args.surprisal_csv} grammar={args.grammar_csv}", flush=True)


if __name__ == "__main__":
    main()
