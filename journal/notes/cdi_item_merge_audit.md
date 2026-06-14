# CDI item-merge audit (io-proc datasets)

Audit of where CDI items are lost across the merges feeding the joint io-proc
model, and a single human-checkable master key. (2026-06-13, task #14.)

## The funnel — and the good news

| step | what happens | items lost |
|---|---|---|
| source CDI per dataset | WS/WG forms: 360–681 items/dataset | — |
| short-code → item_definition | `cdi_short_code_map_ws.csv` (680 codes) / `_wg.csv` (396) | only **non-word columns** (e.g. totlot 804 raw cols → 680 items; the 124 are metadata/headers/dup `feet`) |
| item_definition → CHILDES prob | join to `io_pooled` freq lookup (752-item universe) | **0** — every dataset's items match cleanly |
| **200-item subsample** | stratified class × difficulty for tractability | **552 of 752** (by design) |

**So there is no accidental item loss** — the merges are clean. The only
reduction is the *intentional* 200-item subsample.

## RESOLVED 2026-06-13 — dropped the subsample, now use ALL items

MCF: "aren't we using the full Words & Sentences 680 items? We certainly should
be, not subsampling." Correct. The 200-item subsample is gone:
- `prepare_proc_dp_bundle.R` now defaults to `n_items="all"` (Inf → no subsample;
  pass a number for a quick test fit). proc_dp: **J=681**, N=504,656.
- Joint bundle regenerated: **J=681** (was 200), **N=738,695** (was 217,796, 3.4×),
  I=348, both-channel still 97. SEEDLings now contributes its full WG set (≈360),
  BabyView 604/681 — vs the old 107/165.
- **Compute:** N×3.4 ⇒ the joint fit is ~12–17 h (was ~5 h). Over Sherlock's 16 h
  limit → run on **GCP (uncapped)**, recover posterior offline from CSVs (the
  post-sampling step OOMs at this N, as before).

### On the "752" universe (MCF asked: are these really words?)
752 = **WS(680) ∪ WG-extra(~72)** — a union across forms, all real CDI vocabulary
categories (nouns 352, predicates 175, function_words 107, **other 118**). NOT
grammar/complexity items. The "other" 118 are the CDI's sound-effects
(`vroom`, `grrr`, `cockadoodledoo`), games/routines (`patty cake`,
`this little piggy`), and misc words (`friend`, `please`) — all on the checklist.
**Two debatable items to consider excluding:** `child's own name` (always
produced → uninformative Rasch item) and the multi-word routines. Left in for now
(matches CDI scoring); trivial to drop if we want.

## Artifact

`data/peekbank/cdi_master_item_key.csv` (built by `model/scripts/build_cdi_master_key.R`):
one row per CHILDES-prob item × {short-code, CHILDES prob, lexical class,
in_model_200, n_datasets, presence in each of the 6 datasets}, sorted with the
model items first. This is the single human-checkable merge key MCF asked for —
eyeball the short-code↔item rows and the per-dataset coverage. **72 items lack a
WS short-code** (WG-only or io-derived) — flagged in the `short` column as NA for
review.

## Endgame (per MCF)

A single mapping-key spreadsheet, consistent across datasets — this CSV is the
first version. Next: fold in the WG short-codes for the 72 NAs, and (if we add
Mahr&Edwards / SEEDLings-LWL) extend the per-dataset columns. Then decide the
subsample size as a deliberate, documented choice.
