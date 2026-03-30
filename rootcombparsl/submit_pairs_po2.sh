#!/bin/bash
# Submit Po2 curve jobs in pairs with PBS dependencies so at most 2 geometries
# run simultaneously. Each pair waits for the previous pair to finish.
#
# Usage:
#   bash submit_pairs_po2.sh          # all 25
#   bash submit_pairs_po2.sh 0 9      # geom_000..geom_009 only
#
# Dependencies: requires PBS qsub with -W depend=afterany support.

CURVE_DIR="/dodrio/scratch/projects/2025_060/Joachim/po2_curve"
BATCH_SIZE=2

START=${1:-0}
END=${2:-24}

echo "Submitting geom_$(printf '%03d' ${START})..geom_$(printf '%03d' ${END}) in pairs of ${BATCH_SIZE}"
echo ""

PREV_DEP=""   # dependency string for the next pair

i=${START}
while [ ${i} -le ${END} ]; do
    PAIR_JOBS=""

    for ((k=0; k<BATCH_SIZE && i<=END; k++, i++)); do
        IDX=$(printf "%03d" ${i})
        SCRIPT="${CURVE_DIR}/geom_${IDX}/submit_geom_${IDX}.sh"

        if [ ! -f "${SCRIPT}" ]; then
            echo "  geom_${IDX}: SKIPPED (script not found)"
            continue
        fi

        if [ -n "${PREV_DEP}" ]; then
            JOBID=$(qsub -W depend=afterany:${PREV_DEP} "${SCRIPT}")
        else
            JOBID=$(qsub "${SCRIPT}")
        fi

        JOBID_CLEAN="${JOBID%%.*}"   # strip .pbs... suffix if present
        echo "  geom_${IDX}: ${JOBID}"
        PAIR_JOBS="${PAIR_JOBS}${PAIR_JOBS:+:}${JOBID_CLEAN}"
    done

    PREV_DEP="${PAIR_JOBS}"
done

echo ""
echo "Done. Monitor with: bash check_po2_curve.sh"
echo "      or:           qstat -u vsc45694"
