## Training trajectories (24M rung, 3 seeds)

The `development/` models are each run's best-validation endpoint only. For the largest budget
(24M words) we also release **intermediate checkpoints along training**, retrained with the
identical configuration for **seeds 42, 0 and 123**:

```
trajectory/seed{S}/rung24M/
  step-NNNNNN/          # weights + tokenizer + config + state.json (step, epoch, why, LR)
  log_history.json      # Trainer log: train/val loss per step/epoch, best checkpoint, saved steps
  surprisal.csv         # per-CDI-word mean surprisal (nats) at every saved step
```

- **Checkpoint grid:** step 0 (initialization), 60 log-spaced steps up to 12 epochs (67,464
  steps), and every epoch end, truncated at the early-stopping point: 68–70 checkpoints per
  seed. One epoch is ~5,620 steps (batch 8 × 1024 tokens).
- **Best vs stop:** validation loss bottoms out at epoch 6 (seed 42) or 7 (seeds 0, 123);
  training stops 4 epochs later (step 56,220 / 61,798 / 61,831). **Checkpoints after the best
  epoch are past the validation minimum (overfitting)**; `log_history.json` gives
  `best_model_checkpoint` for each seed.
- **Learning rate:** linear decay is defined over the 40-epoch maximum, so at the stopping
  point the LR is still ~75% of its 1e-4 peak. `state.json` records the LR at each checkpoint.
- **Fidelity:** each seed's best checkpoint reproduces the published 24M ladder model's per-word
  surprisal at r = 0.984–0.988 (mean within 0.04 nats), with best validation loss within
  0.003 nats of the original runs.
- `surprisal.csv` uses the same CDI contexts and scoring as the ladder analysis (128-token left
  context, up to 50 contexts per word; weights stored in bf16, scored in fp32).

```python
sub = "trajectory/seed42/rung24M/step-005622"           # end of epoch 1
model = AutoModelForCausalLM.from_pretrained("mcxfrank/childes-gpt2-ladder", subfolder=sub, torch_dtype="auto")
```

