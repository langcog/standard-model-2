"""Warm-start continuation: load a saved CHILDES checkpoint and continue on the nonce corpus at a
constant low LR, logging per-NONCE surprisal (acquisition) + per-CDI-word surprisal (known-word
control) along a shared log-spaced step axis. No-shuffle => the staggered cohort order is preserved.
"""
import argparse, os, sys
import torch
from datasets import load_dataset
from transformers import (AutoTokenizer, GPT2LMHeadModel, DataCollatorForLanguageModeling,
                          Trainer, TrainingArguments, set_seed)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from surprisal_callback import WordSurprisalCallback
from nonce_surprisal_callback import NonceSurprisalCallback


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init_from", required=True, help="saved checkpoint dir (warm start)")
    ap.add_argument("--train_file", required=True)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--nonce_probe", required=True)
    ap.add_argument("--nonce_csv", required=True)
    ap.add_argument("--cdi_contexts_jsonl", required=True)
    ap.add_argument("--surprisal_csv", required=True)
    ap.add_argument("--num_train_epochs", type=float, default=1.0)
    ap.add_argument("--per_device_batch_size", type=int, default=16)
    ap.add_argument("--learning_rate", type=float, default=1e-5)
    ap.add_argument("--lr_scheduler_type", default="constant")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n_points", type=int, default=60)
    ap.add_argument("--seq_len", type=int, default=1024)
    args = ap.parse_args()
    set_seed(args.seed)

    tok = AutoTokenizer.from_pretrained(args.tokenizer_dir)
    if tok.pad_token is None: tok.pad_token = tok.eos_token
    print(f"[continue] warm-start from {args.init_from}", flush=True)
    model = GPT2LMHeadModel.from_pretrained(args.init_from)
    print(f"[continue] params {sum(p.numel() for p in model.parameters()):,}", flush=True)

    raw = load_dataset("text", data_files={"train": args.train_file}, keep_linebreaks=True)
    tokenized = raw.map(lambda ex: tok(ex["text"]), batched=True, remove_columns=["text"],
                        desc="tok", num_proc=4)
    bs = args.seq_len
    def group(ex):
        cat = {k: sum(ex[k], []) for k in ex}
        n = (len(cat["input_ids"]) // bs) * bs
        r = {k: [t[i:i+bs] for i in range(0, n, bs)] for k, t in cat.items()}
        r["labels"] = r["input_ids"].copy(); return r
    lm = tokenized.map(group, batched=True, num_proc=4, desc="group")
    print(f"[continue] train blocks {len(lm['train'])}", flush=True)

    targs = TrainingArguments(
        output_dir=args.output_dir, overwrite_output_dir=True, do_train=True, do_eval=False,
        num_train_epochs=args.num_train_epochs, per_device_train_batch_size=args.per_device_batch_size,
        learning_rate=args.learning_rate, lr_scheduler_type=args.lr_scheduler_type, warmup_steps=0,
        weight_decay=0.0, adam_beta1=0.9, adam_beta2=0.999, adam_epsilon=1e-8,
        seed=args.seed, data_seed=args.seed, logging_steps=50, eval_strategy="no", save_strategy="no",
        report_to=[], bf16=torch.cuda.is_bf16_supported(), fp16=not torch.cuda.is_bf16_supported(),
        dataloader_num_workers=0)   # kernel 5.4 hangs with worker procs (HF warning); data is small

    class NoShuffle(Trainer):
        def _get_train_sampler(self): return torch.utils.data.SequentialSampler(self.train_dataset)

    cbs = [NonceSurprisalCallback(args.nonce_probe, args.nonce_csv, n_points=args.n_points),
           # known-word control: sparse (8 pts) + small subsample -- just to show known words stay flat
           WordSurprisalCallback(contexts_jsonl=args.cdi_contexts_jsonl, out_csv=args.surprisal_csv,
                                 n_points=8, eval_batch_size=128, max_ctx=128, max_per_word=8)]
    trainer = NoShuffle(model=model, args=targs, train_dataset=lm["train"], tokenizer=tok,
                        data_collator=DataCollatorForLanguageModeling(tokenizer=tok, mlm=False), callbacks=cbs)
    trainer.train()
    print(f"[continue] DONE nonce={args.nonce_csv}", flush=True)


if __name__ == "__main__":
    main()
