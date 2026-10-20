#!/usr/bin/env bash
set -euo pipefail
source "$(cd "$(dirname "$0")" && pwd)/hpc_environment.sh"

echo "Foreground transport run; progress stays on this terminal and is saved under Results_transport_SI/."
python3 -u pfas_runner.py transport "$@" \
  2>&1 | tee -a transport_workflow.log
