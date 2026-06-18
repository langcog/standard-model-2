"""
Per-CDI-word surprisal probe for the autoresearch (nanochat) training loop.

Logs mean NLL (nats) of each CDI target token given a fixed left context, at
log-spaced training steps, to an append-only CSV. This is the Chang & Bergen
per-word acquisition-curve measurement on the TRAINING-STEP axis: fitting a
sigmoid to mean_nll vs log(step) per word gives that word's acquisition slope.

Contexts come from extract_cdi_contexts_ar.py (fixed length K, target = last token).
CSV columns: step,epoch,word,n_occurrences,mean_nll

The model is called as model(ids) -> logits [B,K,V] (nanochat forward with
targets=None); no attention_mask is needed because every context is length K.
"""
import csv, json, math, os
from collections import defaultdict
import torch


def log_spaced_steps(total_steps, n_points=25, min_step=1):
    if total_steps <= 0:
        return []
    n_points = max(2, n_points)
    lo, hi = math.log(max(1, min_step)), math.log(total_steps)
    raw = [math.exp(lo + (hi - lo) * i / (n_points - 1)) for i in range(n_points)]
    return sorted(set(int(round(x)) for x in raw if 1 <= int(round(x)) <= total_steps))


class CdiProbe:
    def __init__(self, contexts_jsonl, out_csv, total_steps, n_points=25, batch_size=512):
        by_word = defaultdict(list)
        K = None
        with open(contexts_jsonl) as f:
            for line in f:
                if not line.strip():
                    continue
                r = json.loads(line)
                ids = r["ids"]
                if K is None:
                    K = len(ids)
                if len(ids) != K:          # enforce uniform length (no padding needed)
                    continue
                by_word[r["word"]].append(ids)
        self.K = K
        self.words = sorted(by_word)
        self.flat_ids, self.flat_w = [], []
        for wi, w in enumerate(self.words):
            for ids in by_word[w]:
                self.flat_ids.append(ids)
                self.flat_w.append(wi)
        self.out_csv = out_csv
        self.batch_size = batch_size
        self.targets = set(log_spaced_steps(total_steps, n_points))
        self.fired = set()
        self._init_csv()
        tg = sorted(self.targets)
        print(f"[cdi_probe] {len(self.words)} words, {len(self.flat_ids)} contexts, K={K}, "
              f"{len(self.targets)} step targets: {tg[:8]}...{tg[-3:]}", flush=True)

    def _init_csv(self):
        d = os.path.dirname(self.out_csv)
        if d:
            os.makedirs(d, exist_ok=True)
        with open(self.out_csv, "w", newline="") as f:
            csv.writer(f).writerow(["step", "epoch", "word", "n_occurrences", "mean_nll"])

    @torch.no_grad()
    def _evaluate(self, model, step, epoch):
        device = next(model.parameters()).device
        # If the model is torch.compiled, call the underlying eager module so the
        # probe's (B,K) shape doesn't trigger a recompile of the training graph.
        fwd = getattr(model, "_orig_mod", model)
        was_training = model.training
        model.eval()
        K, bs = self.K, self.batch_size
        sums, counts = defaultdict(float), defaultdict(int)
        for s in range(0, len(self.flat_ids), bs):
            ids = torch.tensor(self.flat_ids[s:s + bs], dtype=torch.long, device=device)  # [B,K]
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                logits = fwd(ids)                                # [B,K,V]
            lp = torch.log_softmax(logits[:, K - 2, :].float(), dim=-1)  # predict pos K-1
            nll = -lp.gather(1, ids[:, K - 1:K]).squeeze(1)              # [B]
            nll = nll.detach().cpu().tolist()
            for j, wi in enumerate(self.flat_w[s:s + bs]):
                w = self.words[wi]
                sums[w] += nll[j]
                counts[w] += 1
        with open(self.out_csv, "a", newline="") as f:
            wr = csv.writer(f)
            for w in self.words:
                n = counts[w]
                if n:
                    wr.writerow([step, f"{epoch:.4f}", w, n, f"{sums[w] / n:.6f}"])
        if was_training:
            model.train()

    def maybe(self, model, step, epoch):
        if step in self.targets and step not in self.fired:
            self._evaluate(model, step, epoch)
            self.fired.add(step)

    def final(self, model, step, epoch):
        if step not in self.fired:
            self._evaluate(model, step, epoch)
            self.fired.add(step)
