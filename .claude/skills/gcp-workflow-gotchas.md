---
name: gcp-workflow-gotchas
description: Short list of remote-VM workflow footguns that keep biting. Use as a sanity check before kill/relaunch operations, before scp'ing a script to a running VM, or when an SSH session keeps dropping.
---

# GCP / SSH workflow gotchas

## Don't scp a script that's currently executing

`gcloud compute scp gcp/run_fit.sh sm2-fit-01:...` while a `run_fit.sh` process is mid-execution will corrupt the running shell. bash re-reads its source at line boundaries; if the file changed (even just byte length), the next read can land mid-token and throw "unexpected EOF" or "syntax error".

This bricked the EN family driver between sampling and extracts on 2026-05-18. Either:
- Wait for the current invocation to finish, OR
- Make changes to a separately-named file (`gcp/run_fit_v2.sh`) and switch the caller.

## pkill -9 of session ancestors kills the SSH session

If you `gcloud compute ssh ... --command="pkill -9 -f log_irt_long"` and the SSH session itself is somewhere in the process tree being walked, ssh dies. You'll see:
```
ERROR: (gcloud.compute.ssh) [/usr/bin/ssh] exited with return code [255].
```

This is harmless — the pkill ran, just before SSH could send a clean exit. Reconnect to verify the kill worked:
```
gcloud compute ssh <vm> --command="ps -eo cmd | grep <pattern> | grep -v grep || echo clean"
```

## Detached background launches

The robust pattern for "fire and forget" on a remote VM:

```
gcloud compute ssh <vm> --command="cd ~/path && \
  ENV_VAR=val nohup bash run_thing.sh > thing.log 2>&1 < /dev/null & disown"
```

All four pieces matter:

- `nohup` — survives SIGHUP when ssh disconnects
- `> log 2>&1` — captures all output (otherwise SSH waits on it)
- `< /dev/null` — closes stdin (otherwise SSH waits on input)
- `& disown` — backgrounds AND removes from bash's job table

Without `< /dev/null` and `& disown`, the SSH session can still hang or the process can die when ssh closes.

## Background tasks vs. harness notifications

`Bash(..., run_in_background: true)` (the harness feature) is for tracking a command that returns useful output. For pure "fire on remote and forget" launches, use the inline `nohup ... & disown` pattern above and don't background the Bash call itself — the launcher exits in <1s anyway. The harness will then notify when the launcher returns, which is meaningless (it's just confirming launch). Real progress lives in `gcp_*.log` files on the VM.

For things like "watch this fit until it completes" (genuinely long-running), `ScheduleWakeup` with a 15–30 min delay is more useful than a background Bash, because the harness re-invokes you and you can check log + processes fresh.

## The summaries-symlink trap

`gcp/run_fit.sh` sets `SCRATCH=$PWD/_local` and symlinks `$SCRATCH/standard_model_2/summaries` → `$PWD/fits/summaries` so extract scripts (which default to writing under `$SCRATCH`) land in the right place.

An older version of the script did `mkdir -p $SCRATCH/standard_model_2/summaries` before the symlink check, which created a real directory that swallowed extract output. The current version explicitly migrates and replaces if a real dir is found. If you find extract files missing from `fits/summaries/`, check `_local/standard_model_2/summaries/` for the real dir.

## scp with --recurse and multiple sources

```
gcloud compute scp file1.sh file2.sh file3.R vm:dest/ --recurse
```
flattens all sources into `dest/`, breaking the local path hierarchy. To preserve `gcp/run_fit.sh → vm:standard_model_2/gcp/run_fit.sh`, do one scp per file with the explicit destination path. Annoying but reliable.
