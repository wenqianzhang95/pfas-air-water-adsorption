#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/hpc_environment.sh"

requested_case="${1:-all}"
if [[ $# -gt 0 ]]; then shift; fi

if [[ "$requested_case" == "all" ]]; then
  cases=(1 2 3 4)
elif [[ "$requested_case" =~ ^[1-4]$ ]]; then
  cases=("$requested_case")
else
  echo "Usage: bash run_ssl_live.sh [all|1|2|3|4] [runner options]" >&2
  exit 2
fi

for case_id in "${cases[@]}"; do
  echo "Starting inverse SSL case $case_id"
  python3 -u pfas_runner.py ssl --case "$case_id" "$@" \
    2>&1 | tee -a "ssl_case_${case_id}_workflow.log"
done

