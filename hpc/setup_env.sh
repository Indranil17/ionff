#!/bin/bash
# Create the conda environment on a PBS cluster (run once, on a login node or through hpc/setup_env.pbs).
# Usage: STORAGE=/path/to/storage bash hpc/setup_env.sh
# CONDA_OVERRIDE_CUDA makes the solver pick the CUDA build of OpenMM although the login node has no GPU.
set -euo pipefail

: "${STORAGE:?set STORAGE to a directory with enough quota, e.g. your storage home}"
ENV_DIR="${IONFF_ENV:-$STORAGE/miniforge3/envs/ionff}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"

if [ ! -d "$STORAGE/miniforge3" ]; then
    [ -f "$STORAGE/miniforge.sh" ] || curl -L -o "$STORAGE/miniforge.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    bash "$STORAGE/miniforge.sh" -b -p "$STORAGE/miniforge3"
fi
source "$STORAGE/miniforge3/etc/profile.d/conda.sh"

CONDA_OVERRIDE_CUDA="${CONDA_OVERRIDE_CUDA:-12.4}" conda env create -p "$ENV_DIR" -f "$REPO/environment.yml"
conda activate "$ENV_DIR"
pip install -e "$REPO[test]"
python -c "import psi4, openmm; print('psi4', psi4.__version__, 'openmm', openmm.__version__)"
python -m pytest "$REPO/tests" -q -m "not slow"
