#!/bin/bash
# Run each search-direction prompt through agy (Gemini 3.8 Flash, sandboxed) in parallel.
# Usage: bash run_flash.sh <round_dir>
R=$1
for p in "$R"/prompts/d*.txt; do
  n=$(basename "$p" .txt)
  ( cat "$R/prompts/_context.txt" "$p" | { prompt=$(cat); timeout 1500 agy --model gemini-3.8-flash-high --sandbox -p="$prompt"; } > "$R/out/$n.txt" 2> "$R/out/$n.err" ) &
done
wait
