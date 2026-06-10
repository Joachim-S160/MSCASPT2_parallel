#!/bin/bash
# Po2 parallel SO-CASSI — master submission script.
#
# Submits the full 4-phase PBS dependency chain:
#   Phase 1: Setup (GATEWAY+SEWARD+RASSCF×3) — 1 job
#   Phase 2: Workers (CASPT2(Only=N) per root, all spins) — 42 jobs at k=5
#   Phase 3: Assembly (collect_heff.py + CASPT2(EFFE) per spin) — 3 jobs
#   Phase 4: RASSI (SO-CASSI) — 1 job
#
# Usage:
#   bash submit_socassi.sh
#   K=10 bash submit_socassi.sh
#   K=5 KILL_INTRUDERS=0.80 bash submit_socassi.sh
#
# All parameters are env-var overridable for future Python+PBS integration.
# If shared/ already contains all 13 files, setup is skipped.

set -euo pipefail

# ---------------------------------------------------------------------------
# Parameters (overridable via env vars)
# ---------------------------------------------------------------------------
K="${K:-5}"
KILL_INTRUDERS="${KILL_INTRUDERS:-0.70}"

declare -A SPIN_ROOTS=([quintet]=80 [triplet]=90 [singlet]=40)

echo "========================================"
echo "  Po2 parallel SO-CASSI submission"
echo "  K=${K}  KILL_INTRUDERS=${KILL_INTRUDERS}"
echo "========================================"
echo "Date: $(date)"
echo ""

# ---------------------------------------------------------------------------
# Phase 1: Setup (only if shared/ not already populated)
# ---------------------------------------------------------------------------
SHARED_OK=0
if [ -d shared ]; then
    N_IPH=$(ls shared/rasscf_*.JobIph 2>/dev/null | wc -l)
    N_SEW=$(ls shared/seward.* 2>/dev/null | wc -l)
    if [ "${N_IPH}" -ge 3 ] && [ "${N_SEW}" -ge 10 ]; then
        SHARED_OK=1
    fi
fi

if [ "${SHARED_OK}" -eq 1 ]; then
    echo "shared/ already populated ($(ls shared/ | wc -l) files) — skipping setup."
    DEP_SETUP=""
else
    SETUP_ID=$(qsub run_setup.pbs | cut -d. -f1)
    echo "Phase 1 — Setup submitted: ${SETUP_ID}"
    DEP_SETUP="-W depend=afterok:${SETUP_ID}"
fi
echo ""

# ---------------------------------------------------------------------------
# Phase 2 + 3: Workers and assembly per spin
# ---------------------------------------------------------------------------
ALL_ASSEMBLE_IDS=""

for SPIN in quintet triplet singlet; do
    MAX_ROOT="${SPIN_ROOTS[${SPIN}]}"
    BATCH_IDX=0
    WORKER_ID_LIST=""
    N_WORKERS=0

    for START in $(seq 1 "${K}" "${MAX_ROOT}"); do
        END=$(( START + K - 1 ))
        [ "${END}" -gt "${MAX_ROOT}" ] && END="${MAX_ROOT}"

        JOB_ID=$(qsub \
            ${DEP_SETUP} \
            -v "K=${K},SPIN=${SPIN},ROOT_START=${START},ROOT_END=${END},BATCH_IDX=${BATCH_IDX}" \
            run_worker.pbs | cut -d. -f1)

        WORKER_ID_LIST="${WORKER_ID_LIST}:${JOB_ID}"
        BATCH_IDX=$(( BATCH_IDX + 1 ))
        N_WORKERS=$(( N_WORKERS + 1 ))
    done

    # Assembly waits for all workers of this spin
    ASSEMBLE_ID=$(qsub \
        -W "depend=afterok${WORKER_ID_LIST}" \
        -v "SPIN=${SPIN},N_ROOTS=${MAX_ROOT},KILL_INTRUDERS=${KILL_INTRUDERS}" \
        run_assemble.pbs | cut -d. -f1)

    echo "Phase 2/3 — ${SPIN}: ${N_WORKERS} workers → assembly ${ASSEMBLE_ID}"
    ALL_ASSEMBLE_IDS="${ALL_ASSEMBLE_IDS}:${ASSEMBLE_ID}"
done

echo ""

# ---------------------------------------------------------------------------
# Phase 4: RASSI — waits for all 3 assembly jobs
# ---------------------------------------------------------------------------
RASSI_ID=$(qsub \
    -W "depend=afterok${ALL_ASSEMBLE_IDS}" \
    run_rassi.pbs | cut -d. -f1)

echo "Phase 4 — RASSI submitted: ${RASSI_ID}"
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
TOTAL_WORKERS=$(( (80 + K - 1) / K + (90 + K - 1) / K + (40 + K - 1) / K ))
TOTAL_JOBS=$(( 1 + TOTAL_WORKERS + 3 + 1 ))

echo "========================================"
echo "  Submission complete"
echo "  k=${K}  KILL_INTRUDERS=${KILL_INTRUDERS}"
echo "  Workers: ${TOTAL_WORKERS}  Total jobs: ${TOTAL_JOBS}"
echo ""
echo "  Monitor:   qstat -u \$USER"
echo "  After run: python3 collect_heff.py --spin <spin> --n-roots <N>"
echo "             (auto-run by assembly jobs — check heff_*.json)"
echo "             Check rassi.log for SO-CASSI state energies."
echo "========================================"
