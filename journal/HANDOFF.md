# HANDOFF — Standard Model 2 (2026-06-11)

State snapshot to resume cleanly after a long session. Newest state wins; verify
against git/GCP before acting.

## Where master is (origin/master @ bf82beb)
All merged via PRs #23–#37:
- **Repo reorg** (#34): `studies/`, `cluster/`, `figs/`, `reports/`, `journal/`; `outputs/` gone. Map in `/MOVES.md`.
- **build_cache fix** (#35): Fig 4 uses the recent 682-item EN psi; reorg paths repointed.
- **Table 1** (#30): full dataset inventory, landscape page.
- **Fig demographics** (#26, #36): cross-sectional + longitudinal composite.
- **Fig 3 redesign + EN/NO fans** (#37, MERGED): 6 panels — (A) English + (B) Norwegian
  *model-implied* input fans, (C) observed io fan, (D) processing fan on top; (E) factor×channel
  variance partition, (F) σ_r sensitivity below.

## OPEN / STRANDED — pick up here

### 1. Panel-E addition is ORPHANED (re-home it)
Commit **`ee03396`** ("Fig 3 panel E: add EN/NO imputed input-efficiency shares") is on branch
**`fig3-channel-partition`** on origin but **NOT in master** — PR #37 already merged+closed, and a
later push re-created the branch as an orphan ("[new branch]" message was the tell). Work is safe,
just unmerged. It adds the imputed π_α points (io 2.7% / NO 3.9% / EN 7.7%, observed=filled vs
imputed=open) to panel E. To land: new branch off origin/master → apply `ee03396`'s diff
(`build_input_cache.R` partition `source`/`kind` cols + EN/NO rows; `standard_model.qmd` pPart
`shape`/dodge; regen `fig3_input.rds`) → PR. No conflict with MCF's `f9b2496` (that only completed
the io-results prose sentence). Local worktree `../sm2-fig3` holds it.

### 2. NO D fit — likely a GCP refit never pulled
Local NO summary is **May/pre-dedup** vintage (draws 05-23, loo 05-03; psi 673 items vs the
dedup'd **736**-item bundle — same stale-psi gap EN had). MCF believes a NO D refit ran on GCP
~06-10. Both nodes TERMINATED. **ACTION:** start `sm2-fit-02` (project `hs-babyview`, zone
`us-central1-a`; NO is on the 400 GB node), look in `~/standard_model_2/fits/summaries/` for a
recent `long_no_freq_slopes_norwegian_psi.csv` / `.summary` (parallel to the EN psi pulled from
`-01`), `gcloud compute scp` it down. Then the panel-E NO value (provisional 3.9%) + any NO Fig-4
work update. Stop the node after.

### 3. xsec uncap refit = TASK 5 (compute nearly done)
Running via **`nohup`** (logs/xsec_uncapped.log), **30/31** (finishing Cantonese). NOT
harness-tracked → no auto-notification; an exit line `[[ 00_build EXIT=… ]]` lands in the log on
completion. When `fits.rds` reassembles:
1. sanity-check `cross_sectional_demographics/cache/fits.rds` (xsec Ns = true archive counts: EN 8685, NO 7358…; meta k 31/17);
2. rerun `paper/build_table1.R`; drop the "≤1,200" caption clause;
3. compare demographics meta capped-vs-uncapped (expect tighter CIs, same signs);
4. **move `cross_sectional_demographics/` → `studies/`** (deferred reorg piece) + fix the 2 paper paths (`standard_model.qmd` source/readRDS, `build_table1.R` frames);
5. re-render; PR.

## Gotchas (cost real time this session)
- **Don't `pkill 00_build.R` before relaunching** — it kills the healthy running job. Most "deaths" were self-inflicted + the harness `run_in_background` mechanism (Falcon/lid were NOT it). `nohup`+`disown` survives.
- **Don't push to a branch whose PR already merged** — GitHub auto-deletes it on merge; a push re-creates an orphan. Check `gh pr view` state first.
- **`build_input_cache.R` needs `fits/`** (gitignored) → run from the MAIN repo, not a worktree; copy the regenerated `fig3_input.rds` into the worktree, then `git checkout HEAD` the main copy.
- **Stale psi**: EN/NO/proc fits had May psi vs dedup'd bundles. Always check psi item-count == bundle word_info count before using `delta_j`.

## Other active worktrees/sessions (possible parallel edits)
`../sm2-fig3` (orphan panel-E), `../standard_model_2-demometa` [fig-demographics-inline-meta],
`.claude/worktrees/disjoint-analysis` [claude/disjoint-control]. Coordinate before editing the qmd.
