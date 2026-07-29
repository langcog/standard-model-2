#!/bin/bash
D=/data2/mcfrank/babylm
mkdir -p "$D"; cd "$D" || exit 1
LOG="$D/.stage.log"; rm -f "$D/.stage.DONE"
echo "START $(date)" > "$LOG"
curl -sL "https://osf.io/download/661518e7219e71402df6a816/" -o train_100M.zip 2>> "$LOG"
echo "zip size: $(ls -lh train_100M.zip | awk '{print $5}')" >> "$LOG"
unzip -o train_100M.zip -d . >> "$LOG" 2>&1
echo "== per-source word counts ==" >> "$LOG"
find . -name "*.train" -o -name "*.txt" | while read f; do echo "$(wc -w < "$f") $f" >> "$LOG"; done
echo "DONE $(date)" >> "$LOG"
echo done > "$D/.stage.DONE"
