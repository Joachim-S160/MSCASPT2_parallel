#!/bin/bash
# Po2 intruder-state test — master submission script.
#
# Submits the full 5-phase PBS DAG for all N_GEOM geometries:
#
#   Phase 1:  N_GEOM setup jobs (parallel, independent)
#   Phase 2:  N_GEOM × (ceil(Q/K)+ceil(T/K)+ceil(S/K)) worker jobs
#             (each geometry's workers depend on that geometry's setup)
#   Phase 3:  1 global consensus job (depends on ALL workers)
#   Phase 4:  N_GEOM × 3 assembly jobs (all depend on consensus)
#   Phase 5:  N_GEOM RASSI jobs (each depends on its 3 assembly jobs)
#
# Usage:
#   cd MSCASPT2_parallel/SO_CASSI_test/po2_intruder_test/
#   bash submit_intruder_test.sh [--dry-run]
#
# --dry-run: print all qsub commands without submitting.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    echo "=== DRY RUN MODE (no jobs submitted) ==="
fi

qsub_or_dry() {
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "DRYRUN: qsub $*" >&2
        # Return a fake job ID for dependency tracking
        echo "FAKE.$(( RANDOM ))"
    else
        qsub "$@" | tr -d '\n'
    fi
}

echo "========================================================"
echo "  Po2 intruder-state test — DAG submission"
echo "  N_GEOM=${N_GEOM}  K=${K}  KILL_INTRUDERS=${KILL_INTRUDERS}"
echo "  CIROOT: Q=${CIROOT_Q}  T=${CIROOT_T}  S=${CIROOT_S}"
echo "  Script dir: ${SCRIPT_DIR}"
echo "========================================================"
echo ""

# ---------------------------------------------------------------------------
# Phase 1: Setup (one job per geometry)
# ---------------------------------------------------------------------------
echo "--- Phase 1: Setup ---"
SETUP_IDS=()
for i in $(seq 0 $((N_GEOM - 1))); do
    GEOM_PADDED=$(printf "%02d" "${i}")
    XYZ_NAME="${GEOM_LIST[$i]}"
    XYZ_FILE="${SCRIPT_DIR}/xyz/${XYZ_NAME}.xyz"

    if [ ! -f "${XYZ_FILE}" ]; then
        echo "ERROR: xyz file not found: ${XYZ_FILE}"
        exit 1
    fi

    mkdir -p "${SCRIPT_DIR}/run_${GEOM_PADDED}/shared"

    JID=$(qsub_or_dry \
        -N "po2_intruder_setup_g${i}" \
        -v "GEOM_IDX=${i},XYZ_FILE=${XYZ_FILE}" \
        "${SCRIPT_DIR}/run_setup.pbs")
    SETUP_IDS+=("${JID}")
    echo "  g${i} (${XYZ_NAME}): setup → ${JID}"
done
echo ""

# ---------------------------------------------------------------------------
# Phase 2: Workers (per geometry × per spin batch)
# ---------------------------------------------------------------------------
echo "--- Phase 2: Workers ---"
ALL_WORKER_IDS=()
declare -A SPIN_ROOTS=([quintet]="${CIROOT_Q}" [triplet]="${CIROOT_T}" [singlet]="${CIROOT_S}")

for i in $(seq 0 $((N_GEOM - 1))); do
    SETUP_DEP="-W depend=afterok:${SETUP_IDS[$i]}"
    GEOM_WORKER_IDS=()

    for SPIN in quintet triplet singlet; do
        MAX_ROOT="${SPIN_ROOTS[$SPIN]}"
        BATCH_IDX=0

        for START in $(seq 1 "${K}" "${MAX_ROOT}"); do
            END=$(( START + K - 1 ))
            [ "${END}" -gt "${MAX_ROOT}" ] && END="${MAX_ROOT}"

            JID=$(qsub_or_dry ${SETUP_DEP} \
                -N "po2_intruder_worker_g${i}_${SPIN:0:1}${BATCH_IDX}" \
                -v "GEOM_IDX=${i},K=${K},SPIN=${SPIN},ROOT_START=${START},ROOT_END=${END},BATCH_IDX=${BATCH_IDX},FROZEN=${FROZEN}" \
                "${SCRIPT_DIR}/run_worker.pbs")

            GEOM_WORKER_IDS+=("${JID}")
            ALL_WORKER_IDS+=("${JID}")
            BATCH_IDX=$(( BATCH_IDX + 1 ))
        done
    done

    N_WORKERS_THIS="${#GEOM_WORKER_IDS[@]}"
    echo "  g${i}: ${N_WORKERS_THIS} workers submitted (depend on ${SETUP_IDS[$i]})"
done

echo "  Total workers: ${#ALL_WORKER_IDS[@]}"
echo ""

# ---------------------------------------------------------------------------
# Phase 3: Global consensus (depends on ALL workers)
# ---------------------------------------------------------------------------
echo "--- Phase 3: Global consensus ---"
WORKER_DEP_STR="afterok:$(IFS=':'; echo "${ALL_WORKER_IDS[*]}")"

CONSENSUS_ID=$(qsub_or_dry \
    -N "po2_intruder_consensus" \
    -W "depend=${WORKER_DEP_STR}" \
    "${SCRIPT_DIR}/run_global_consensus.pbs")
echo "  Consensus job: ${CONSENSUS_ID}  (depends on ${#ALL_WORKER_IDS[@]} workers)"
echo ""

# ---------------------------------------------------------------------------
# Phase 4: Assembly (per geometry × 3 spins, all depend on consensus)
# ---------------------------------------------------------------------------
echo "--- Phase 4: Assembly ---"
declare -A ASSEMBLE_IDS_BY_GEOM
CONSENSUS_DEP="-W depend=afterok:${CONSENSUS_ID}"

for i in $(seq 0 $((N_GEOM - 1))); do
    GEOM_ASSEMBLE_IDS=()
    declare -A SPIN_N_ROOTS=([quintet]="${CIROOT_Q}" [triplet]="${CIROOT_T}" [singlet]="${CIROOT_S}")

    for SPIN in quintet triplet singlet; do
        N_ROOTS="${SPIN_N_ROOTS[$SPIN]}"
        JID=$(qsub_or_dry ${CONSENSUS_DEP} \
            -N "po2_intruder_assemble_g${i}_${SPIN:0:1}" \
            -v "GEOM_IDX=${i},SPIN=${SPIN},N_ROOTS=${N_ROOTS}" \
            "${SCRIPT_DIR}/run_assemble.pbs")
        GEOM_ASSEMBLE_IDS+=("${JID}")
    done

    ASSEMBLE_IDS_BY_GEOM[$i]="${GEOM_ASSEMBLE_IDS[*]}"
    echo "  g${i}: assembly jobs: ${GEOM_ASSEMBLE_IDS[*]}"
done
echo ""

# ---------------------------------------------------------------------------
# Phase 5: RASSI (per geometry, depends on its 3 assembly jobs)
# ---------------------------------------------------------------------------
echo "--- Phase 5: RASSI ---"
for i in $(seq 0 $((N_GEOM - 1))); do
    ASSEMBLE_DEP_STR="afterok:$(echo "${ASSEMBLE_IDS_BY_GEOM[$i]}" | tr ' ' ':')"

    JID=$(qsub_or_dry \
        -N "po2_intruder_rassi_g${i}" \
        -W "depend=${ASSEMBLE_DEP_STR}" \
        -v "GEOM_IDX=${i}" \
        "${SCRIPT_DIR}/run_rassi.pbs")
    echo "  g${i}: RASSI → ${JID}  (depends on ${ASSEMBLE_IDS_BY_GEOM[$i]})"
done
echo ""

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
N_Q_BATCHES=$(( (CIROOT_Q + K - 1) / K ))
N_T_BATCHES=$(( (CIROOT_T + K - 1) / K ))
N_S_BATCHES=$(( (CIROOT_S + K - 1) / K ))
N_WORKERS_PER_GEOM=$(( N_Q_BATCHES + N_T_BATCHES + N_S_BATCHES ))
N_TOTAL_WORKERS=$(( N_WORKERS_PER_GEOM * N_GEOM ))

echo "========================================================"
echo "  Job count summary"
echo "  Setup:     ${N_GEOM}"
echo "  Workers:   ${N_TOTAL_WORKERS}  (${N_WORKERS_PER_GEOM}/geom: ${N_Q_BATCHES}Q+${N_T_BATCHES}T+${N_S_BATCHES}S)"
echo "  Consensus: 1"
echo "  Assembly:  $(( N_GEOM * 3 ))  (3/geom)"
echo "  RASSI:     ${N_GEOM}"
echo "  TOTAL:     $(( N_GEOM + N_TOTAL_WORKERS + 1 + N_GEOM * 3 + N_GEOM ))"
echo "========================================================"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  (dry run — no jobs were submitted)"
fi
