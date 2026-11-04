# classic_chunking — does a chunking model explain v2's word difficulty?  ·  exploratory

**Self-contained study; not (yet) in the paper.** Question: Rowland et al. (2026, ch. 4, Language 0–5)
explain vocabulary growth with CLASSIC, a chunking model of the *input*. Does CLASSIC's per-word
learning time predict v2's per-word difficulty δ_j **beyond frequency**? (v2's RQ2 null: frequency adds
nothing beyond a free δ_j.)

## Result (main run: 40 streams × 60k CHILDES caregiver utterances, ages 10–36 mo)

10-fold CV R² for predicting item difficulty; `log_t` = mean log utterance index at which CLASSIC first
stores the item as a whole-word chunk; class = v2 lexical class (4 levels).

| Model | δ_j (named fit, n = 195) | empirical difficulty (n = 676) |
|---|---|---|
| class only | 0.18 | 0.16 |
| class × log_p (class-specific slopes) | 0.21 | 0.25 |
| class × log_p + log_p² (best frequency baseline) | 0.40 | 0.37 |
| class × log_t (CLASSIC alone) | 0.43 | 0.42 |
| best frequency baseline + log_t | **0.52** | **0.46** |

Partial r(log_t, difficulty | class, log_p) = +0.62 [0.51, 0.71] (δ_j), +0.53 [0.46, 0.59] (empirical).
Tokens of the item heard before learning (`log_tok`) adds almost nothing (partial r ≤ 0.09), so the signal
is CLASSIC's learning time, not raw exposure.

Sensitivity (CV R² of best frequency baseline + log_t; baselines are 0.398 / 0.373 throughout):

| Run | δ_j | empirical |
|---|---|---|
| `main` 60k utts, ages 10–36 | 0.521 | 0.461 |
| `n15k` 15k utts | 0.494 | 0.419 |
| `n120k` 120k utts | 0.485 | 0.460 |
| `age60` ages 10–60 | 0.495 | 0.462 |
| `lr20` learning rate 0.2 | 0.497 | 0.466 |

CLASSIC alone is 0.39–0.43 (δ_j) and 0.37–0.42 (empirical): comparable to the best frequency baseline,
not clearly better. The robust claim is the *increment* (+0.05 to +0.12), not CLASSIC replacing frequency.

## What carries the increment? (`06_features.py`)

Most of CLASSIC's gain over frequency is reproduced by cheap distributional features of the same
CHILDES input. Per item, from the 10–36 mo pool: share of occurrences that are utterance-final
(`p_final`), share that are the whole utterance (`p_alone`), log mean length of the utterances containing
it (`log_len`), and entropy of its left/right neighbours (`ctx_ent`); plus phoneme count (`n_ph`).
CV R², main run (B = class × log_p + log_p²):

| Model | δ_j (n = 195) | empirical (n = 676) |
|---|---|---|
| B | 0.398 | 0.373 |
| B + CLASSIC `log_t` | 0.521 | 0.461 |
| B + features | 0.485 | 0.484 |
| B + features + `log_t` | 0.532 | 0.522 |
| class × `log_t` + features, no frequency | 0.490 | 0.506 |

- The four features add +0.087 (δ_j) and +0.111 (empirical) to B, against +0.123 and +0.088 for `log_t`:
  about 70% of CLASSIC's increment on δ_j and all of it on the empirical target. A plain feature model
  without CLASSIC matches CLASSIC alone.
- `log_len` is the one that matters: words that occur in longer utterances are harder (Spearman with δ_j
  +0.49; dropping it costs −0.02 to −0.03 CV R², dropping any other feature changes R² by < 0.01).
- What is left for CLASSIC after features + `n_ph`: +0.049 / +0.033 (δ_j / empirical, main) and
  +0.03 to +0.04 in the other settings, except `n15k` empirical (+0.008). Partial
  r(`log_t`, difficulty | class, log_p, features, n_ph) = +0.52 [0.41, 0.63] (δ_j), +0.46 [0.39, 0.52]
  (empirical), so a residual exists but is small and is not shown to be chunking-specific.
- Frequency measure: CHILDES token counts (`log_cnt`) behave like v2's `prob` as the baseline (0.37 / 0.41).
- Caveats: four features chosen after seeing the problem, no feature × class interactions, n = 195 for
  δ_j, and the features are from the same corpus CLASSIC was run on.

Reading: v2's word difficulty depends on utterance context (length, position, diversity) in a way a
frequency count misses; CLASSIC is one way of encoding that, not obviously a better one.

Calibration (main run, median CDI words known): 89 @5% of the stream, 195 @15%, 392 @50%, 498 @100%.
Rowland et al. report 104 / 373 / 558 at 5 / 13 / 48% for British children: same ballpark, not a replication.

## Pipeline

`sh studies/classic_chunking/run_all.sh /path/to/cmudict.dict` (from the repo root; ~1 min after the
one-off CHILDES pull). Intermediate data go to `work/` (gitignored, ~76 MB); small outputs to `results/`.

| Step | Script | Output |
|---|---|---|
| v2 item difficulty table | `00_export_v2_items.R` | `work/v2_items_delta.csv` |
| CHILDES Eng-NA caregiver utterances (childesr) | `01_pull_childes.R` | `work/childes_eng_na_utts.csv.gz` |
| Phonemize utterances + CDI items (CMUdict) | `02_prep_data.py` | `work/prep.pkl` |
| CLASSIC runs (`classic_fast.py`) | `03_run_experiment.py` | `work/runs_<tag>.csv`, `growth_<tag>.csv` |
| Correlations, partials, class-FE R², bootstrap CIs | `analyze.py` | `results/<tag>_correlations.txt` |
| Cross-validated R² for baselines vs CLASSIC | `05_cv_baselines.py` | `results/<tag>_cv_r2.txt` |
| Distributional features vs the CLASSIC increment | `06_features.py` (+ `cvlib.py`) | `results/<tag>_features.txt` |

`results/item_table_main.csv` has the per-item table (δ_j, empirical difficulty, log_p, log_t, …) for the main run.

## Provenance and assumptions

- **CLASSIC** (Jones et al.; Rowland et al. 2026). `classic_fast.py` re-implements the OSF script
  `classic.py` (https://osf.io/9s4d7, file https://osf.io/download/wvc6y/) rule for rule with set/dict
  lookups. Details not stated in the chapter but taken from the script: access probability
  1/(1+exp(−(4.5−x−1))), learning rate 0.10 per adjacent accessed pair, parse order
  multiword → lexical → sublexical → phonemic, a word is "known" once its string is a lexical or
  multiword chunk. The OSF script was read, not executed, and is not vendored here.
- **Input is American, not British**: CHILDES `Eng-NA`, speaker ≠ target child, ages 10–60 mo, matching
  the CDI norms behind v2 rather than the Language 0–5 corpus. Streams are age-ordered samples.
- **Phonemes**: CMUdict ARPAbet, stress stripped, first pronunciation only, a few hand-added entries and
  two-word compound splitting for OOV tokens. 91.1% of utterances survive (names, Spanish, codes dropped;
  the chapter keeps 96%). 678/681 CDI items have phonemes (child/pet/babysitter names dropped).
- **Difficulty targets.** `delta198` is v2's per-word ψ/δ_j from `figs/longitudinal/psi_freq_regression_per_word.csv`
  (197 items; 195 heard in the streams). That file's `jj` is a private 1..198 index, so the join is **by item
  name**. `emp_diff` = −logit(production rate, ages 22–28), r ≈ 0.98 with δ_j where both exist.
  **Do not use `fits/summaries/long_no_freq_slopes_psi.csv` for item-level comparisons**: its 682 rows do
  not align with the 681-item bundle (correlation with production ≈ 0).

## Not done

- Sub-lexical/phonological structure and feature × class interactions as sources of CLASSIC's residual
  increment (position, length and context are covered above).
- Child level: whether CLASSIC can produce v2's efficiency/acceleration structure, and CLASSIC run at the
  empirical input spread (σ_r ≈ 0.4) against v2's 4–7% input share.
- Only 195 items have a valid named δ_j; the all-item test uses empirical difficulty. Full 681-item δ_j needs the
  stale `long_no_freq_slopes_psi.csv` regenerated: see [`journal/notes/TODO_psi_681_refit.md`](../../journal/notes/TODO_psi_681_refit.md).
