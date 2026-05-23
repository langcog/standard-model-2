# standard_model_2 project skills

Skills here are auto-discoverable by Claude Code when you're working in
this repo. Their frontmatter `description` fields are matched against
your request to decide when to load them.

## What's here

- **`gcp-stan-fitting.md`** — Infrastructure playbook: VM config, threading
  rule, compile flags, disk requirements, env vars, entry points. Read
  this before launching anything.
- **`stan-mixing-failures.md`** — Diagnose & recover from non-convergent
  fits. Covers the warm-start / boundary-parking pattern that hit us on
  EN demo_alpha (Rhat 2.7, n_eff 5), and why bumping adapt_delta /
  max_treedepth is not the answer.
- **`psis-loo-large-datasets.md`** — Chunked PSIS-LOO recipe for log_lik
  matrices >~30 GB. The vanilla `loo()` segfaults at this scale;
  `sherlock/extract_loo_thinned.R` ships the fix.
- **`gcp-workflow-gotchas.md`** — Short list of remote-VM footguns:
  don't scp into running scripts, SSH-dies-on-pkill, detached-launch
  pattern, summaries-symlink trap.

## How to invoke

The Skill tool: `Skill(skill="gcp-stan-fitting")`. Or just describe what
you're doing — the description field matches request intent. Adding a
new skill: create another `.md` file in this directory with the same
frontmatter format.

## Relationship to global skills

Generic versions of these lessons live in `~/.claude/skills/`:
`gcp-large-stan`, `sherlock-stan-fitting`, `stan-mixing-cold-start`,
`psis-loo-chunked`, `remote-vm-gotchas`. Those generalize across
projects; the ones in this directory have `standard_model_2`-specific
VM names (`sm2-fit-01/02`), repo paths, and variant lists baked in.

When both a local and global skill exist with the same name (e.g.,
`gcp-stan-fitting`), local takes precedence within this repo.
