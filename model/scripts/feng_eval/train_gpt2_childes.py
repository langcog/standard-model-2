"""
Train GPT-2 small from scratch on CHILDES (Feng et al. 2026 specification),
logging per-CDI-word mean surprisal at log-spaced training steps.

Matches Feng et al. 2026 §B settings for the CHILDES condition (GPT-2 small,
LR 1e-4 linear, no warmup, 20 epochs, batch 8/GPU, no in-epoch shuffling,
Adam β=(0.9, 0.999), ε=1e-8, seq length 1024). The data file
`CHILDES_train_ordered.txt` is consumed line-by-line in the given order.

Eval set for the callback is the precomputed CDI-context JSONL.

Usage:
    python train_gpt2_childes.py \
        --train_file CHILDES_train_ordered.txt \
        --val_file CHILDES_val_ordered.txt \
        --tokenizer_dir tokenizers/GPT2_CHILDES \
        --config_file tokenizers/GPT2-small_config/config.json \
        --output_dir /tmp/run_seed42 \
        --cdi_contexts_jsonl cdi_contexts.jsonl \
        --surprisal_csv surprisal_seed42.csv \
        --num_train_epochs 20 \
        --per_device_batch_size 8 \
        --learning_rate 1e-4 \
        --seed 42 \
        --n_log_points 80
"""

import argparse
import os
import sys

import torch
from datasets import load_dataset
from transformers import (
    AutoTokenizer,
    GPT2Config,
    GPT2LMHeadModel,
    DataCollatorForLanguageModeling,
    Trainer,
    TrainingArguments,
    set_seed,
)

# Allow running from anywhere; ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from surprisal_callback import WordSurprisalCallback


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_file", required=True)
    ap.add_argument("--val_file", required=True)
    ap.add_argument("--tokenizer_dir", required=True)
    ap.add_argument("--config_file", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--cdi_contexts_jsonl", required=True)
    ap.add_argument("--surprisal_csv", required=True)
    ap.add_argument("--num_train_epochs", type=float, default=20.0)
    ap.add_argument("--per_device_batch_size", type=int, default=8)
    ap.add_argument("--learning_rate", type=float, default=1e-4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n_log_points", type=int, default=80)
    ap.add_argument("--save_total_limit", type=int, default=1)
    ap.add_argument("--seq_len", type=int, default=1024)
    ap.add_argument("--logging_steps", type=int, default=50)
    ap.add_argument("--eval_callback_batch_size", type=int, default=64)
    ap.add_argument("--no_save", action="store_true",
                    help="Skip saving any model checkpoints (we don't need them)")
    args = ap.parse_args()

    set_seed(args.seed)

    print(f"[train] loading tokenizer from {args.tokenizer_dir}", flush=True)
    tokenizer = AutoTokenizer.from_pretrained(args.tokenizer_dir)
    # GPT-2 has no pad token by default — use eos for padding in collator.
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    print(f"[train] loading config from {args.config_file}", flush=True)
    config = GPT2Config.from_pretrained(args.config_file)
    # Make config's vocab match the tokenizer
    config.vocab_size = tokenizer.vocab_size
    config.bos_token_id = tokenizer.bos_token_id if tokenizer.bos_token_id is not None else tokenizer.eos_token_id
    config.eos_token_id = tokenizer.eos_token_id

    print(f"[train] init GPT-2 from config (vocab {config.vocab_size}, "
          f"n_layer {config.n_layer}, n_embd {config.n_embd})", flush=True)
    model = GPT2LMHeadModel(config)
    n_params = sum(p.numel() for p in model.parameters())
    print(f"[train] model params: {n_params:,}", flush=True)

    # ---- Datasets ----
    print(f"[train] loading text datasets", flush=True)
    raw = load_dataset(
        "text",
        data_files={"train": args.train_file, "validation": args.val_file},
        keep_linebreaks=True,
    )
    print(f"[train] raw train lines: {len(raw['train'])}, val lines: {len(raw['validation'])}",
          flush=True)

    def tokenize_function(examples):
        return tokenizer(examples["text"])

    tokenized = raw.map(
        tokenize_function,
        batched=True,
        remove_columns=["text"],
        desc="Tokenizing",
        num_proc=4,
    )

    block_size = args.seq_len

    def group_texts(examples):
        # Concatenate all texts; chunk into blocks of block_size.
        # Mirrors HF run_clm.py.
        concatenated = {k: sum(examples[k], []) for k in examples.keys()}
        total_length = len(concatenated[list(examples.keys())[0]])
        total_length = (total_length // block_size) * block_size
        result = {
            k: [t[i:i + block_size] for i in range(0, total_length, block_size)]
            for k, t in concatenated.items()
        }
        result["labels"] = result["input_ids"].copy()
        return result

    lm_dataset = tokenized.map(
        group_texts,
        batched=True,
        num_proc=4,
        desc=f"Grouping into {block_size}-token blocks",
    )
    print(f"[train] grouped train blocks: {len(lm_dataset['train'])}, "
          f"val blocks: {len(lm_dataset['validation'])}", flush=True)

    # ---- Training ----
    save_strategy = "no" if args.no_save else "epoch"
    targs = TrainingArguments(
        output_dir=args.output_dir,
        overwrite_output_dir=True,
        do_train=True,
        do_eval=True,
        num_train_epochs=args.num_train_epochs,
        per_device_train_batch_size=args.per_device_batch_size,
        per_device_eval_batch_size=args.per_device_batch_size,
        learning_rate=args.learning_rate,
        lr_scheduler_type="linear",
        warmup_steps=0,
        weight_decay=0.0,
        adam_beta1=0.9,
        adam_beta2=0.999,
        adam_epsilon=1e-8,
        seed=args.seed,
        data_seed=args.seed,
        logging_steps=args.logging_steps,
        evaluation_strategy="epoch",
        save_strategy=save_strategy,
        save_total_limit=args.save_total_limit,
        # match Feng et al.: no shuffling within epoch
        # HF Trainer arg name is `dataloader_drop_last`/`dataloader_shuffle`;
        # the shuffling flag controls in-epoch shuffling for train loader.
        report_to=[],  # no W&B
        bf16=torch.cuda.is_bf16_supported(),
        fp16=not torch.cuda.is_bf16_supported(),
        dataloader_num_workers=2,
        ddp_find_unused_parameters=False,
    )
    # `train_dataloader_shuffle` was added in transformers≥4.41. For 4.39 we
    # subclass Trainer to override get_train_dataloader and disable shuffling.

    class NoShuffleTrainer(Trainer):
        def _get_train_sampler(self):
            # Sequential sampler => no shuffling
            return torch.utils.data.SequentialSampler(self.train_dataset)

    callback = WordSurprisalCallback(
        contexts_jsonl=args.cdi_contexts_jsonl,
        out_csv=args.surprisal_csv,
        n_points=args.n_log_points,
        eval_batch_size=args.eval_callback_batch_size,
        max_ctx=block_size,
    )

    trainer = NoShuffleTrainer(
        model=model,
        args=targs,
        train_dataset=lm_dataset["train"],
        eval_dataset=lm_dataset["validation"],
        tokenizer=tokenizer,
        data_collator=DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False),
        callbacks=[callback],
    )

    trainer.train()
    if not args.no_save:
        trainer.save_model(os.path.join(args.output_dir, "final"))
    print(f"[train] DONE. surprisal log: {args.surprisal_csv}", flush=True)


if __name__ == "__main__":
    main()
