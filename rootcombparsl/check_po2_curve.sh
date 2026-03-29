#!/bin/bash
# Check progress of the 25-geometry Po2 dissociation curve run.

CURVE_DIR="/dodrio/scratch/projects/2025_060/Joachim/po2_curve"

done_count=0
running_count=0
pending_count=0

for i in $(seq 0 24); do
    IDX=$(printf "%03d" $i)
    RASSI_LOG="${CURVE_DIR}/geom_${IDX}/final_rassi/final_rassi.log"
    ANY_LOG=$(ls "${CURVE_DIR}/geom_${IDX}"/singlet/rasscf/rasscf_full.log 2>/dev/null || true)

    if grep -q "Happy landing" "${RASSI_LOG}" 2>/dev/null; then
        STATUS="DONE"
        ((done_count++))
    elif [ -n "${ANY_LOG}" ] || ls "${CURVE_DIR}/geom_${IDX}"/singlet/root*/root*.log 2>/dev/null | head -1 | grep -q .; then
        STATUS="RUNNING"
        ((running_count++))
    else
        STATUS="pending"
        ((pending_count++))
    fi

    printf "  geom_%s  %s\n" "${IDX}" "${STATUS}"
done

echo ""
echo "Summary: ${done_count}/25 done, ${running_count} running, ${pending_count} pending"
