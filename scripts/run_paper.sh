#!/usr/bin/env bash
# Reproduce all DFT calculations shown in Gül & Roth, Physica A 674, 130686 (2025).
#
#   scripts/run_paper.sh [outdir]          (default: data/runs)
#
# Environment: JOBS (parallel runs, default 4), OMP_NUM_THREADS (threads per run, default 2).
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:-$root/data/runs}"
bin="$root/src/sw_interface"
jobs="${JOBS:-4}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

[ -x "$bin" ] || { echo "build the solver first: make -C src" >&2; exit 1; }
mkdir -p "$out"

# "q eps2" per line; "sym" solves for the eps2 with rho2^l = rho2^v.
cases=(
  "1.0 0.10" "1.0 0.30" "1.0 sym" "1.0 0.60" "1.0 1.20"   # Figs. 1 and 2
  "1.5 0.60" "1.5 sym" "1.5 1.00"                          # Fig. 3(a)
  "2.0 1.20" "2.0 sym" "2.0 1.50"                          # Fig. 3(b), Fig. 4
  "2.5 1.80" "2.5 sym" "2.5 2.20"                          # Fig. 3(c)
  "3.0 2.00" "3.0 sym" "3.0 3.00"                          # Fig. 3(d)
)

run_case() {
  local q=$1 e=$2 log
  log="$out/log_q${q}_${e}.txt"
  if [ "$e" = sym ]; then
    "$bin" -q "$q" -s -e 1.0 -o "$out" >"$log" 2>&1
  else
    "$bin" -q "$q" -e "$e" -o "$out" >"$log" 2>&1
  fi
  echo "done: q=$q eps2=$e  $(grep gamma "$log")"
}
export -f run_case
export bin out

printf '%s\n' "${cases[@]}" | xargs -P "$jobs" -L 1 bash -c 'run_case $0 $1'
