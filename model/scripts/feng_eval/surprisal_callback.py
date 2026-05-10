"""
HuggingFace Trainer callback: log per-CDI-word mean surprisal at log-spaced
training steps.

Loaded by run_clm_with_surprisal_logging.py (a thin wrapper around
run_clm_no_shuffling.py).

Output CSV (appended per logged step):
    step,word,n_occurrences,mean_nll,sum_nll
where mean_nll is the average -log p(target | left context) over all stored
occurrences of `word`, in NATS (natural log).
"""

import csv
import json
import os
from collections import defaultdict
from typing import List

import torch
from transformers import TrainerCallback


def log_spaced_steps(total_steps: int, n_points: int = 80, min_step: int = 1) -> List[int]:
    """Log-spaced integer step targets in [min_step, total_steps], inclusive."""
    import math
    if total_steps <= 0:
        return []
    n_points = max(2, n_points)
    log_lo = math.log(max(1, min_step))
    log_hi = math.log(total_steps)
    raw = [math.exp(log_lo + (log_hi - log_lo) * i / (n_points - 1))
           for i in range(n_points)]
    return sorted(set(int(round(x)) for x in raw if 1 <= int(round(x)) <= total_steps))


class WordSurprisalCallback(TrainerCallback):
    """At each scheduled training step, compute per-CDI-word mean NLL on a
    precomputed eval set and append rows to a CSV.

    Parameters
    ----------
    contexts_jsonl : path to JSONL produced by extract_cdi_contexts.py
    out_csv : output CSV path (appended)
    n_points : number of log-spaced steps to log at
    eval_batch_size : how many contexts to score per forward pass
    max_ctx : truncate context to last `max_ctx` tokens
    """

    def __init__(self, contexts_jsonl, out_csv, n_points=80,
                 eval_batch_size=64, max_ctx=1024, also_log_step_0=True):
        self.contexts_jsonl = contexts_jsonl
        self.out_csv = out_csv
        self.n_points = n_points
        self.eval_batch_size = eval_batch_size
        self.max_ctx = max_ctx
        self.also_log_step_0 = also_log_step_0

        # Load contexts grouped by word
        self.by_word = defaultdict(list)  # word -> list[(ctx_ids, target_pos)]
        with open(contexts_jsonl) as f:
            for line in f:
                if not line.strip():
                    continue
                r = json.loads(line)
                ctx = r["ctx"]
                if len(ctx) > self.max_ctx:
                    ctx = ctx[-self.max_ctx:]
                # target was at last position of original ctx; after truncation
                # it's still at last position.
                self.by_word[r["word"]].append((ctx, len(ctx) - 1))

        self.words = sorted(self.by_word.keys())
        self.n_contexts = sum(len(v) for v in self.by_word.values())

        # Flat list of (word_idx, ctx_ids, tgt_pos) for batching
        self._flat = []
        for wi, w in enumerate(self.words):
            for ctx, pos in self.by_word[w]:
                self._flat.append((wi, ctx, pos))

        self._step_targets = None
        self._fired = set()
        self._csv_initialized = False

    def _init_csv(self):
        os.makedirs(os.path.dirname(self.out_csv) or ".", exist_ok=True)
        if not os.path.exists(self.out_csv):
            with open(self.out_csv, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["step", "epoch", "word", "n_occurrences",
                            "mean_nll", "sum_nll"])
        self._csv_initialized = True

    def _maybe_set_targets(self, max_steps):
        if self._step_targets is None:
            self._step_targets = set(log_spaced_steps(max_steps, self.n_points))
            print(f"[surprisal_callback] log-spaced step targets ({len(self._step_targets)}): "
                  f"{sorted(self._step_targets)[:10]} ... "
                  f"{sorted(self._step_targets)[-5:]}", flush=True)

    @torch.no_grad()
    def _evaluate(self, model, device, step, epoch):
        model.eval()
        # For each context, score the target. We do this by feeding the prefix
        # ctx[:pos+1] (with the target as the LAST token) and reading the
        # log-prob of ctx[pos] given ctx[:pos]. That requires a forward pass
        # where input_ids = ctx, and we look at logits at position pos-1 for the
        # token at position pos (standard causal LM scoring).
        sums = defaultdict(float)
        counts = defaultdict(int)

        bs = self.eval_batch_size
        for start in range(0, len(self._flat), bs):
            batch = self._flat[start:start + bs]
            max_len = max(len(c) for _, c, _ in batch)
            ids = torch.full((len(batch), max_len),
                             fill_value=0, dtype=torch.long, device=device)
            attn = torch.zeros((len(batch), max_len), dtype=torch.long, device=device)
            for i, (_, ctx, _) in enumerate(batch):
                L = len(ctx)
                ids[i, :L] = torch.tensor(ctx, device=device)
                attn[i, :L] = 1

            out = model(input_ids=ids, attention_mask=attn, use_cache=False)
            logits = out.logits  # [B, T, V]
            log_probs = torch.log_softmax(logits.float(), dim=-1)

            for i, (wi, ctx, pos) in enumerate(batch):
                # nll of token at position `pos` given tokens [0, pos-1].
                # Logits at index pos-1 predict the token at index pos.
                # Skip pos == 0 (no left context).
                if pos < 1:
                    continue
                tgt = ctx[pos]
                lp = log_probs[i, pos - 1, tgt].item()
                w = self.words[wi]
                sums[w] += -lp
                counts[w] += 1

        # Append rows
        with open(self.out_csv, "a", newline="") as f:
            wr = csv.writer(f)
            for w in self.words:
                n = counts[w]
                if n == 0:
                    continue
                mean_nll = sums[w] / n
                wr.writerow([step, epoch, w, n, f"{mean_nll:.6f}", f"{sums[w]:.6f}"])
        model.train()

    def on_train_begin(self, args, state, control, **kwargs):
        if not self._csv_initialized:
            self._init_csv()
        # max_steps may not yet be finalized at train begin; we re-set on step.

    def on_step_end(self, args, state, control, **kwargs):
        self._maybe_set_targets(state.max_steps)
        s = state.global_step
        if s in self._fired:
            return
        if s in self._step_targets or (self.also_log_step_0 and s == 1 and 1 not in self._fired):
            # On first step, also log the (near-)pre-training state by using the
            # current model (which is essentially un-trained after step 1).
            model = kwargs.get("model")
            if model is None:
                return
            device = next(model.parameters()).device
            self._evaluate(model, device, step=s, epoch=state.epoch)
            self._fired.add(s)

    def on_train_end(self, args, state, control, **kwargs):
        # Final step, in case it wasn't in the schedule
        s = state.global_step
        if s in self._fired:
            return
        model = kwargs.get("model")
        if model is None:
            return
        device = next(model.parameters()).device
        self._evaluate(model, device, step=s, epoch=state.epoch)
        self._fired.add(s)
