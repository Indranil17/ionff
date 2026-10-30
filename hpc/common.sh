# Shared setup sourced by every job script. Set STORAGE (and optionally IONFF_ENV, IONFF_REPO) before qsub,
# for example with qsub -v STORAGE=...,IONFF_REPO=...
: "${STORAGE:?STORAGE not set}"
ENV_DIR="${IONFF_ENV:-$STORAGE/miniforge3/envs/ionff}"
REPO="${IONFF_REPO:-${PBS_O_WORKDIR:-$PWD}}"
source "$STORAGE/miniforge3/etc/profile.d/conda.sh"
conda activate "$ENV_DIR"
cd "$REPO"

# psi4 scratch on node-local disk when the scheduler provides it
export PSI_SCRATCH="${SCRATCHDIR:-${TMPDIR:-/tmp}}/psi4_$$"
mkdir -p "$PSI_SCRATCH"
# PBS presets NCPUS=1, so take the requested count from PBS_NCPUS
export NCPUS="${PBS_NCPUS:-8}"
trap 'rm -rf "$PSI_SCRATCH"' EXIT
