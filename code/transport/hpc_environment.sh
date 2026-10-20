#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PACKAGE_DIR"

if ! command -v matlab >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    module purge
    module load matlab
  fi
fi

if ! command -v matlab >/dev/null 2>&1; then
  echo "ERROR: MATLAB is unavailable. Load a MATLAB module or set PFAS_MATLAB_CMD." >&2
  exit 2
fi
if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is unavailable." >&2
  exit 2
fi

job_tmp="${SLURM_JOBTMP:-${TMPDIR:-/tmp}}"
if [[ -z "${MATLAB_PREFDIR:-}" ]]; then
  MATLAB_PREFDIR="$(mktemp -d "$job_tmp/pfas-matlab-${SLURM_JOB_ID:-local}-XXXXXX")"
  export MATLAB_PREFDIR
fi
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export PYTHONUNBUFFERED=1

