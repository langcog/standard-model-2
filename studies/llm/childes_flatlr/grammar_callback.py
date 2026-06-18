"""HuggingFace Trainer callback: at log-spaced training steps, log grammar accuracy
(Zorro + optionally BLiMP) AND val bits-per-token, so the grammar/competence trajectory
aligns step-for-step with the CDI-surprisal trajectory (WordSurprisalCallback).

Output CSV (appended): step,epoch,val_bpb,zorro,blimp
Grammar eval is the expensive part, so this fires at fewer points (n_points ~10-12) and
subsamples pairs per paradigm (max_pairs). A full eval can be run post-hoc with eval_grammar.py.
"""
import csv, math, os
import torch
from transformers import TrainerCallback

from surprisal_callback import log_spaced_steps   # same dir
from eval_grammar import grammar_accuracy


class GrammarTrajectoryCallback(TrainerCallback):
    def __init__(self, tokenizer, zorro_dir, out_csv, n_points=12, max_pairs=100,
                 do_blimp=True, blimp_configs=None, eval_batch_size=64, max_len=128,
                 val_input_ids=None, also_log_step_1=True):
        self.tok = tokenizer
        self.zorro_dir = zorro_dir
        self.out_csv = out_csv
        self.n_points = n_points
        self.max_pairs = max_pairs
        self.do_blimp = do_blimp
        self.blimp_configs = blimp_configs
        self.bs = eval_batch_size
        self.max_len = max_len
        self.val_input_ids = val_input_ids       # tensor [N, L] subsample, or None
        self.also_log_step_1 = also_log_step_1
        self._targets = None
        self._fired = set()
        self._init = False

    def _init_csv(self):
        os.makedirs(os.path.dirname(self.out_csv) or ".", exist_ok=True)
        with open(self.out_csv, "w", newline="") as f:
            csv.writer(f).writerow(["step", "epoch", "val_bpb", "zorro", "blimp"])
        self._init = True

    @torch.no_grad()
    def _val_bpb(self, model, device):
        # Manual token NLL (NOT model(labels=) -- HF ForCausalLMLoss blows up memory
        # casting/fusing full-sequence logits; the manual gather path is what the
        # grammar scorer uses and is memory-safe).
        if self.val_input_ids is None:
            return ""
        import torch.nn.functional as F
        model.eval()
        tot_nll, tot_tok = 0.0, 0
        ids_all = self.val_input_ids
        vbs = 4
        for i in range(0, ids_all.shape[0], vbs):
            b = ids_all[i:i + vbs].to(device)
            logits = model(input_ids=b).logits                # [B,L,V]
            logp = F.log_softmax(logits.float(), dim=-1)
            tgt = b[:, 1:]
            lp = logp[:, :-1, :].gather(-1, tgt.unsqueeze(-1)).squeeze(-1)   # [B,L-1]
            tot_nll += float((-lp).sum())
            tot_tok += tgt.numel()
        return tot_nll / tot_tok / math.log(2)     # bits per token

    @torch.no_grad()
    def _eval(self, model, device, step, epoch):
        val_bpb = self._val_bpb(model, device)
        gbs = min(self.bs, 32)   # minimal-pair sents <=128 tok: [32,128,V] log_softmax ~safe mid-training
        res = grammar_accuracy(model, self.tok, self.zorro_dir, device,
                               batch_size=gbs, max_pairs=self.max_pairs,
                               max_len=self.max_len, do_blimp=self.do_blimp,
                               blimp_configs=self.blimp_configs)
        blimp = res["blimp_overall"] if res["blimp_overall"] is not None else ""
        with open(self.out_csv, "a", newline="") as f:
            csv.writer(f).writerow([step, f"{epoch:.4f}",
                                    f"{val_bpb:.5f}" if val_bpb != "" else "",
                                    f"{res['zorro_overall']:.5f}",
                                    f"{blimp:.5f}" if blimp != "" else ""])
        print(f"[grammar_traj] step {step} ep {epoch:.2f}: val_bpb={val_bpb} "
              f"zorro={res['zorro_overall']:.3f} blimp={blimp}", flush=True)
        model.train()

    def on_train_begin(self, args, state, control, **kwargs):
        if not self._init:
            self._init_csv()

    def on_step_end(self, args, state, control, **kwargs):
        if self._targets is None and state.max_steps > 0:
            self._targets = set(log_spaced_steps(state.max_steps, self.n_points))
            print(f"[grammar_traj] targets: {sorted(self._targets)}", flush=True)
        s = state.global_step
        if s in self._fired or self._targets is None:
            return
        if s in self._targets or (self.also_log_step_1 and s == 1):
            model = kwargs.get("model")
            if model is None:
                return
            self._eval(model, next(model.parameters()).device, s, state.epoch)
            self._fired.add(s)

    def on_train_end(self, args, state, control, **kwargs):
        s = state.global_step
        if s in self._fired:
            return
        model = kwargs.get("model")
        if model is None:
            return
        self._eval(model, next(model.parameters()).device, s, state.epoch)
        self._fired.add(s)
