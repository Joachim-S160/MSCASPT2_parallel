#!/bin/bash
# Po2 k-scaling benchmark — master submission script.
#
# Usage:
#   bash submit_scaling.sh           # sample mode: 3 batches per k (default)
#   bash submit_scaling.sh --full    # full mode: all 90 triplet roots per k
#
# Submits:
#   1 setup job (run_setup.pbs)
#   k-sweep worker jobs (run_worker.pbs), each depending on the setup job
#
# Sample mode (default): 3 representative batches per k from triplet roots.
#   Root offsets: 1, 11, 21  (span low / mid / upper range of the triplet block)
#   Total: 1 setup + 3 × 4 k-values = 13 jobs
#
# Full mode (--full): all ceil(90/k) batches for k ∈ {1,2,5,10}
#   Total: 1 setup + 90+45+18+9 = 163 jobs
#
# After all jobs complete:
#   python3 collect_scaling.py --plot scaling_results.png

set -euo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPIN="triplet"
MAX_ROOT=90
K_VALUES=(1 2 5 10)
FULL_MODE=0

for arg in "$@"; do
    case "${arg}" in
        --full) FULL_MODE=1 ;;
        *) echo "Unknown argument: ${arg}"; exit 1 ;;
    esac
done

cd "${BENCH_DIR}"

# ---------------------------------------------------------------------------
# Verify shared/ exists or submit setup job
# ---------------------------------------------------------------------------
if [ -d "${BENCH_DIR}/shared" ] && [ -f "${BENCH_DIR}/shared/rasscf_triplet.JobIph" ]; then
    echo "shared/ already exists — skipping setup job."
    echo "To force a fresh setup: rm -rf shared/ and re-run this script."
    SETUP_ID=""
    DEPEND_FLAG=""
else
    echo "Submitting setup job..."
    SETUP_SUBMIT=$(qsub "${BENCH_DIR}/run_setup.pbs")
    SETUP_ID="${SETUP_SUBMIT%%.*}"
    DEPEND_FLAG="-W depend=afterok:${SETUP_ID}"
    echo "  Setup job: ${SETUP_ID}"
fi

# ---------------------------------------------------------------------------
# Submit worker jobs
# ---------------------------------------------------------------------------
TOTAL_WORKERS=0

for K in "${K_VALUES[@]}"; do
    BATCH_IDX=0

    if [ "${FULL_MODE}" -eq 1 ]; then
        # Full sweep: all roots 1..MAX_ROOT in steps of K
        ROOT=1
        while [ "${ROOT}" -le "${MAX_ROOT}" ]; do
            END=$(( ROOT + K - 1 ))
            [ "${END}" -gt "${MAX_ROOT}" ] && END="${MAX_ROOT}"

            JOB_ID=$(qsub \
                ${DEPEND_FLAG} \
                -v "K=${K},SPIN=${SPIN},ROOT_START=${ROOT},ROOT_END=${END},BATCH_IDX=${BATCH_IDX}" \
                "${BENCH_DIR}/run_worker.pbs" \
            )
            echo "  k=${K} b=${BATCH_IDX} roots=${ROOT}-${END} → ${JOB_ID%%.*}"
            BATCH_IDX=$(( BATCH_IDX + 1 ))
            ROOT=$(( END + 1 ))
            TOTAL_WORKERS=$(( TOTAL_WORKERS + 1 ))
        done
    else
        # Sample mode: 3 batches at root offsets 1, 11, 21
        for OFFSET in 1 11 21; do
            END=$(( OFFSET + K - 1 ))
            [ "${END}" -gt "${MAX_ROOT}" ] && END="${MAX_ROOT}"

            JOB_ID=$(qsub \
                ${DEPEND_FLAG} \
                -v "K=${K},SPIN=${SPIN},ROOT_START=${OFFSET},ROOT_END=${END},BATCH_IDX=${BATCH_IDX}" \
                "${BENCH_DIR}/run_worker.pbs" \
            )
            echo "  k=${K} b=${BATCH_IDX} roots=${OFFSET}-${END} → ${JOB_ID%%.*}"
            BATCH_IDX=$(( BATCH_IDX + 1 ))
            TOTAL_WORKERS=$(( TOTAL_WORKERS + 1 ))
        done
    fi
done

echo ""
if [ -n "${SETUP_ID:-}" ]; then
    echo "Submitted: 1 setup (${SETUP_ID}) + ${TOTAL_WORKERS} workers"
    echo "All workers depend on afterok:${SETUP_ID}"
else
    echo "Submitted: ${TOTAL_WORKERS} workers (no setup — shared/ already present)"
fi
echo ""
echo "Monitor: qstat -u \$USER"
echo "After all jobs finish:"
echo "  python3 collect_scaling.py --plot scaling_results.png"
