"""Callback: at log-spaced steps, log summed-subword NLL of each NONCE word across its held-out
probe contexts (acquisition curve). Output: step,epoch,nonce,cohort,intro_frac,train_freq,mean_span_nll,n_ctx
Mirrors surprisal_callback but targets multi-subword nonce spans (sum NLL over the span)."""
import csv, json, os
from collections import defaultdict
import torch, torch.nn.functional as F
from transformers import TrainerCallback
from surprisal_callback import log_spaced_steps


class NonceSurprisalCallback(TrainerCallback):
    def __init__(self, probe_jsonl, out_csv, n_points=60, eval_batch_size=32, also_log_step_1=True):
        self.rows = [json.loads(l) for l in open(probe_jsonl) if l.strip()]
        self.out_csv = out_csv; self.n_points = n_points; self.bs = eval_batch_size
        self.also1 = also_log_step_1
        self.meta = {r["nonce"]: (r["cohort"], r["intro_frac"], r["train_freq"]) for r in self.rows}
        self._t = None; self._fired = set(); self._init = False

    def _init_csv(self):
        os.makedirs(os.path.dirname(self.out_csv) or ".", exist_ok=True)
        with open(self.out_csv, "w", newline="") as f:
            csv.writer(f).writerow(["step","epoch","nonce","cohort","intro_frac","train_freq","mean_span_nll","n_ctx"])
        self._init = True

    @torch.no_grad()
    def _eval(self, model, device, step, epoch):
        model.eval()
        sums, cnts = defaultdict(float), defaultdict(int)
        rows = sorted(self.rows, key=lambda r: len(r["ids"]))
        for i in range(0, len(rows), self.bs):
            batch = rows[i:i + self.bs]; ml = max(len(r["ids"]) for r in batch)
            ids = torch.zeros(len(batch), ml, dtype=torch.long, device=device)
            att = torch.zeros(len(batch), ml, dtype=torch.long, device=device)
            for j, r in enumerate(batch):
                L = len(r["ids"]); ids[j, :L] = torch.tensor(r["ids"], device=device); att[j, :L] = 1
            logp = F.log_softmax(model(input_ids=ids, attention_mask=att).logits.float(), dim=-1)
            # vectorized gather of all span-token logprobs (one sync per batch, no per-token sync)
            rr, pp, tt, owner = [], [], [], []
            for j, r in enumerate(batch):
                s, e = r["span"]
                for t in range(max(s, 1), e):              # token t scored from position t-1
                    rr.append(j); pp.append(t - 1); tt.append(r["ids"][t]); owner.append(j)
            lp = logp[torch.tensor(rr, device=device), torch.tensor(pp, device=device),
                      torch.tensor(tt, device=device)].cpu().tolist()
            for k, j in enumerate(owner):
                sums[batch[j]["nonce"]] += -lp[k]
            for r in batch:
                cnts[r["nonce"]] += 1
        with open(self.out_csv, "a", newline="") as f:
            w = csv.writer(f)
            for nz in sorted(sums):
                c, fr, tf = self.meta[nz]
                w.writerow([step, f"{epoch:.4f}", nz, c, fr, tf, f"{sums[nz]/cnts[nz]:.5f}", cnts[nz]])
        model.train()

    def on_train_begin(self, args, state, control, **kw):
        if not self._init: self._init_csv()

    def on_step_end(self, args, state, control, **kw):
        if self._t is None and state.max_steps > 0:
            self._t = set(log_spaced_steps(state.max_steps, self.n_points))
        s = state.global_step
        if self._t is None or s in self._fired: return
        if s in self._t or (self.also1 and s == 1):
            m = kw.get("model")
            if m is None: return
            self._eval(m, next(m.parameters()).device, s, state.epoch); self._fired.add(s)

    def on_train_end(self, args, state, control, **kw):
        s = state.global_step
        if s in self._fired: return
        m = kw.get("model")
        if m is not None:
            self._eval(m, next(m.parameters()).device, s, state.epoch); self._fired.add(s)
