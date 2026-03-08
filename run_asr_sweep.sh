#!/bin/bash
#SBATCH --job-name=ASR_sweep
#SBATCH --partition=standard
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=80:00:00
#SBATCH --output=/home/%u/AdaptiveASR/logs/job_%x_%j.out
#SBATCH --error=/home/%u/AdaptiveASR/logs/job_%x_%j.err

# ============================================================
# ASR Parameter Sweep — PARAM Kamrupa (IIT Guwahati)
#
# Node spec: 2x Intel Xeon Platinum 8268 = 48 cores, 192 GB RAM
# We request 48 cores and 180 GB (leaving 12 GB headroom for OS).
#
# Parallelism: parfor over all (subject x combo) pairs jointly.
# With 10 subjects x N combos, up to 48 workers run simultaneously.
# Tasks queue behind the 48 active workers — no cores ever idle.
#
# One file per (subject, combo) is written to results/.
# After the sweep completes, accumulateResults merges them
# into one per-subject file (see bottom of launcher).
# ============================================================

set -euo pipefail

PROJECT_ROOT="${HOME}/AdaptiveASR"
LOG_DIR="${PROJECT_ROOT}/logs"
LAUNCHER="run_parameterSweep_launcher"

# ── Sanity checks ─────────────────────────────────────────────
if [[ ! -d "${PROJECT_ROOT}" ]]; then
    echo "[ERROR] Project root not found: ${PROJECT_ROOT}"
    echo "        Edit PROJECT_ROOT in this script to match your cluster path."
    exit 1
fi

mkdir -p "${LOG_DIR}"
cd "${PROJECT_ROOT}"

# ── Environment ───────────────────────────────────────────────
module purge
module load matlab/R2022b

# Pass the allocated CPU count to MATLAB so the launcher opens
# a parpool of exactly the right size without hardcoding it.
export SLURM_WORKERS="${SLURM_CPUS_PER_TASK}"

echo "======================================================"
echo "  Job:     ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"
echo "  Node:    $(hostname)"
echo "  CPUs:    ${SLURM_CPUS_PER_TASK}"
echo "  Memory:  ${SLURM_MEM_PER_NODE} MB"
echo "  Start:   $(date)"
echo "  Root:    ${PROJECT_ROOT}"
echo "======================================================"

# ── Run MATLAB ────────────────────────────────────────────────
matlab -nodesktop -nosplash -r "
try;
    run('${LAUNCHER}.m');
catch e;
    fprintf('[MATLAB FATAL] %s\n', getReport(e, 'extended'));
    exit(1);
end;
exit(0);
"

EXIT_CODE=$?
echo "======================================================"
echo "  End:        $(date)"
echo "  Exit code:  ${EXIT_CODE}"
echo "======================================================"
exit ${EXIT_CODE}
