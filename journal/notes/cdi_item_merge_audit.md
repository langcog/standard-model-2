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

## But two things worth acting on

1. **The subsample leaves 552 items on the table.** 752 CHILDES-prob items exist;
   the model uses 200. More items → better-estimated per-child efficiency (ξ).
   The subsample was a tractability choice; the joint fit is ~5 h, so **raising it
   (400, or all 752) is feasible** and is the cleanest way to recover signal
   ("losing items hurts us"). → candidate robustness/precision lever.
2. **Cross-dataset coverage is uneven.** The 200 chosen items are fully present in
   the four WS/Stanford RT datasets (200/200), but **BabyView 165/200** and
   **SEEDLingS 107/200** (SEEDLingS is the younger WG form). So SEEDLings kids'
   ξ is estimated on only ~107 items. Choosing the subsample to also cover the
   WG/SEEDLings vocabulary (or just enlarging it) would sharpen the input-only
   datasets' efficiency estimates. 341 items are in **all 6** datasets — a natural
   shared core.

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
