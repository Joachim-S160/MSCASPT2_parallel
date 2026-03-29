#!/bin/bash
# Submit Po2 dissociation curve geometry jobs to PBS.
# Usage:
#   bash submit_all_po2.sh          # submit all 25 (geom_000..geom_024)
#   bash submit_all_po2.sh 0 9     # submit geom_000..geom_009 only
#   bash submit_all_po2.sh 10 24   # submit geom_010..geom_024 only

CURVE_DIR="/dodrio/scratch/projects/2025_060/Joachim/po2_curve"

START=${1:-0}
END=${2:-24}

echo "Submitting geom_$(printf '%03d' ${START}) through geom_$(printf '%03d' ${END})..."
echo ""

for i in $(seq ${START} ${END}); do
    IDX=$(printf "%03d" $i)
    SCRIPT="${CURVE_DIR}/geom_${IDX}/submit_geom_${IDX}.sh"
    if [ ! -f "${SCRIPT}" ]; then
        echo "  geom_${IDX}: SKIPPED (script not found — run setup_po2_curve.sh first)"
        continue
    fi
    JOBID=$(qsub "${SCRIPT}")
    echo "  geom_${IDX}: ${JOBID}"
done

echo ""
echo "Done. Monitor with: qstat -u vsc45694"
