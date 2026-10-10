"""
Save bf16 weight snapshots at listed global steps and at every epoch end, into
{output_dir}/trajectory/step-NNNNNN/. Writes outside checkpoint-*, so the Trainer's
save_total_limit rotation and early-stopping best-model logic are untouched.
Weights only (no optimizer state). See ckpt24m/SPEC.md.
"""

import json
import os

import torch
from transformers import TrainerCallback


class TrajectoryCheckpointCallback(TrainerCallback):
    def __init__(self, steps_file, output_dir, tokenizer):
        with open(steps_file) as f:
            self.steps = {int(x) for x in f.read().split()}
        self.root = os.path.join(output_dir, "trajectory")
        self.tokenizer = tokenizer
        self.saved = set()

    def _save(self, model, state, why):
        step = state.global_step
        if step in self.saved:
            return
        d = os.path.join(self.root, f"step-{step:06d}")
        # bf16 copy; the live (training) model is not modified.
        sd = {k: v.detach().to(torch.bfloat16).cpu() for k, v in model.state_dict().items()}
        if getattr(model.config, "tie_word_embeddings", False):
            sd.pop("lm_head.weight", None)  # tied to transformer.wte; restored on load
        model.save_pretrained(d, state_dict=sd, safe_serialization=True)
        cfg = os.path.join(d, "config.json")  # config says fp32 (live model); weights are bf16
        c = json.load(open(cfg))
        c["torch_dtype"] = "bfloat16"
        json.dump(c, open(cfg, "w"), indent=2)
        self.tokenizer.save_pretrained(d)
        lr = next((h["learning_rate"] for h in reversed(state.log_history)
                   if "learning_rate" in h), None)
        with open(os.path.join(d, "state.json"), "w") as f:
            json.dump({"step": step, "epoch": state.epoch, "why": why,
                       "last_logged_lr": lr}, f)
        self.saved.add(step)

    def on_train_begin(self, args, state, control, model=None, **kw):
        if 0 in self.steps:
            self._save(model, state, "init")

    def on_step_end(self, args, state, control, model=None, **kw):
        if state.global_step in self.steps:
            self._save(model, state, "grid")

    def on_epoch_end(self, args, state, control, model=None, **kw):
        self._save(model, state, "epoch_end")

    def on_train_end(self, args, state, control, **kw):
        os.makedirs(self.root, exist_ok=True)
        with open(os.path.join(self.root, "log_history.json"), "w") as f:
            json.dump({"best_model_checkpoint": state.best_model_checkpoint,
                       "best_metric": state.best_metric,
                       "saved_steps": sorted(self.saved),
                       "log_history": state.log_history}, f, indent=1)
