# BabyView raw data

Stanford BabyView corpus: head-mounted body-cam video of children's
daily lives, with auto-transcribed speech and parent-report CDI data.

## Provenance

- Release: BabyView **2025.1**.
- Source: Stanford BabyView project (langcog / BabyView lab).
- Citation: Long, B. L., Sparks, R. Z., Xiang, V., Stojanov, S., Yin, Z.,
  Keene, G., Tan, A. W. M., Feng, S. Y., Nag, A., Zhuang, C.,
  Marchman, V. A., Yamins, D. L. K., & Frank, M. C. (2025). The
  BabyView dataset: High-resolution egocentric videos of infants'
  and young children's everyday experiences. In *Proceedings of the
  Cognitive Computational Neuroscience Conference* (8-page track).

The BabyView release contains:
- Per-video metadata (subject, age, recording length, English content
  fraction).
- Auto-transcribed speech with per-token speaker tags (FEM/MAL/CHI),
  confidence scores, and tokenization.
- Parent-report CDI item-level responses for the longitudinal
  subjects.

## Files

| file | size | what it contains |
|---|---|---|
| `video_metadata_processed.csv` | 1.5 MB | per-video metadata (subject_id, video_id, duration_sec, percent_english, age, ...) |
| `merged_transcripts_parsed.csv` | **835 MB, gitignored** | per-token transcript output (speaker, score, spacy parse) for all videos. Re-derive locally if needed. |
| `cdi_data_oct_2025/babyview-english-{wg,ws}_items.csv` | ~50 KB each | item-level CDI responses (long format) for English admins |
| `cdi_data_oct_2025/babyview-spanish-{wg,ws}_items.csv` | ~50 KB each | item-level CDI responses for Spanish admins (currently unused; pipeline filters to English-dominant subjects) |

## Pipeline

`model/scripts/prepare_babyview.R` consumes:
- video metadata + transcript → per-video adult-token rate
  (`log_r_obs`)
- CDI item-level → admin-level production observations (`y_ij`)

Output: `model/fits/babyview_subset_data.rds` (Stan-ready bundle).
The 835 MB transcript is **not** committed; the prepared bundle is
~500 KB and is committed so Sherlock can fit without re-running prep.

## What's NOT here

- The raw video files (Databrary-restricted).
- Older / newer BabyView releases.
