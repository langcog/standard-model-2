# Research menu #2 — other learning-dynamics & architecture levers for J

Generated 2026-06-18 (overnight) by a background research agent (primary-source verified).
Complements `research_menu.md` (grokking/data). Excludes already-tried: repetition, WD,
grokfast, decay-vs-flat LR.

## Central tension (the artifact control)
Almost every mechanism that curves the loss upward *also* can defer early progress. Three
ways to get J>0: (1) **genuine sequential/compositional emergence** — a new capability
switches on mid-training and compounds (induction heads, saddle-to-saddle rank steps,
"Domino effect"); (2) **schedule-induced reveal** — WSD parks progress in a high-LR
plateau then reveals it at decay (real progress, schedule-imposed timing); (3) **pure
slow-start** — init/warmup/regularization that just delays early gain (the artifact).
**Decisive control = matched-early-token surprisal:** a genuine arm matches/beats baseline
early AND pulls ahead late. Our `returns_score.R` reports `r_early` — J>=0 is REAL only if
r_early stays healthy (~>=0.4) AND final_surp near baseline (~<=6.1). (This is exactly how
grokfast's high-lamb J>0 was caught as an artifact: r_early collapsed to ~0 / negative.)

## Ranked shortlist (G = likely genuine, R = reveal/timing, A = high artifact risk)
| # | lever | train.py knob | genuine? | cost | cite |
|---|---|---|---|---|---|
| 1 | **Depth-over-width** at fixed params | ↑n_layer, ↓n_embd (↓ASPECT_RATIO for deep) | **G** induction/ICL emergence sharper+later w/ depth | mod | [13][1] |
| 2 | **WSD schedule** | warmup→long constant high-LR→short sharp decay | **R→G** verified reveal, not stall | triv | [2] |
| 3 | **Small-init / rich regime** | scale down init (esp. readout) | **G** saddle-to-saddle stepwise head formation | triv-mod | [8][10] |
| 4 | **Curriculum-of-capabilities** | order data foundational→dependent | **G** Domino compounding skill chains | mod (data) | [21][20] |
| 5 | **Large-LR / edge-of-stability** (+QK-norm, already present) | ↑peak LR | **G-lean** catapult into flatter minima/feature learning | triv | [6][11][5] |
| 6 | **μP** | width-scaled init+LR (pairs w/ Muon) | **G enabler** maximal feature learning + LR transfer | mod | [14][10] |
| 7 | **L2-to-preserve-ICL** | small L2 to retain emergent ICL | **G protector** stops ICL transience that kills J | triv | [9] |
| 8 | Sophia optimizer | swap MuonAdamW→Sophia | R/efficiency | mod | [7] |
| 9 | QK-norm | (already in this model) | enabler | triv | [3] |
| 10 | Ramped critical-batch | small batch early→larger late | mild R | mod | [17] |

**Most likely GENUINE:** #1 depth, #3 small-init, #4 curriculum, #5 large-LR — all switch on
a compounding capability rather than deferring early learning. Highest-EV composite:
`depth + small-init + QK-norm(have) + WSD`. Guardrail: #7 L2-preserve-ICL.
Also: fit per-word curves with a BNSL inflection form [22] to tell real slope-change from
monotone-curved stall; re-measure late gains to check ICL transience [9].

## Sources (key)
1 Olsson+ 2022 induction heads https://transformer-circuits.pub/2022/in-context-learning-and-induction-heads/index.html ·
2 Wen+ 2024 WSD river-valley https://arxiv.org/abs/2410.05192 ·
3 QK-norm https://arxiv.org/abs/2511.21377 ·
4 Kosson+ 2024 why warmup https://arxiv.org/abs/2406.09405 ·
5 EoS minimalist https://arxiv.org/abs/2503.02809 ·
6 catapult https://arxiv.org/abs/2003.02218 ·
7 Sophia https://arxiv.org/abs/2305.14342 ·
8 Jacot+ saddle-to-saddle small-init https://arxiv.org/abs/2106.15933 ·
9 Singh+ 2023 ICL transience https://arxiv.org/abs/2311.08360 ·
10 lazy-vs-rich tutorial https://arxiv.org/abs/2404.19719 ·
12 Muon scalable https://arxiv.org/abs/2502.16982 ·
13 Mehta+Gupta 2025 scaling+ICL (depth allocation) https://arxiv.org/abs/2511.06232 ·
14 u-μP https://arxiv.org/html/2407.17465v3 ·
17 McCandlish+ 2018 critical batch https://arxiv.org/abs/1812.06162 ·
20 curriculum learning dynamics https://arxiv.org/abs/2601.21698 ·
21 Liu+ 2025 Physics of Skill Learning (Domino) https://arxiv.org/abs/2501.12391 ·
22 Caballero+ 2023 Broken Neural Scaling Laws https://arxiv.org/abs/2210.14891
