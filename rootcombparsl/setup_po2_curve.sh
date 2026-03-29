#!/bin/bash
# Creates 25 per-geometry run directories under CURVE_DIR.
# Each directory gets config_geom_NNN.yml and submit_geom_NNN.sh.
#
# Pairing: XYZ ascending (po2_000→po2_024), RasOrb descending (system_24→system_0)
# because autoCAS scanned the curve in the reverse direction of the xyz numbering.
#
# Run once on HPC after: git pull (autoCAS4HE) + git pull (MSCASPT2_parallel)

set -e

WORKFLOW_DIR="/dodrio/scratch/projects/starting_2025_097/MSCASPT2_parallel/rootcombparsl"
CURVE_DIR="/dodrio/scratch/projects/2025_060/Joachim/po2_curve"
AUTOCAS_BUILD="/dodrio/scratch/projects/starting_2025_097/autoCAS4HE_built"
XYZ_DIR="${AUTOCAS_BUILD}/autoCAS4HE/tests/molcas/SOCASSI/VQZP_25"
RASORB_DIR="${AUTOCAS_BUILD}/autoCAS4HE/tests/autocas/external_scf/Po2_scaling_test_final_casscf/VQZP/results_N25_rev"

# Verify source paths exist
for D in "${WORKFLOW_DIR}" "${XYZ_DIR}" "${RASORB_DIR}"; do
    if [ ! -d "${D}" ]; then
        echo "ERROR: directory not found: ${D}"
        exit 1
    fi
done

mapfile -t XYZ_FILES < <(ls "${XYZ_DIR}"/po2_*.xyz | sort -V)
mapfile -t RASORB_FILES < <(ls "${RASORB_DIR}"/system_*.RasOrb | sort -rV)

N_XYZ=${#XYZ_FILES[@]}
N_RAS=${#RASORB_FILES[@]}
if [ "${N_XYZ}" -ne 25 ] || [ "${N_RAS}" -ne 25 ]; then
    echo "ERROR: expected 25 xyz and 25 RasOrb files, got ${N_XYZ} xyz and ${N_RAS} RasOrb"
    exit 1
fi

mkdir -p "${CURVE_DIR}"
echo "Setting up ${N_XYZ} geometry directories in ${CURVE_DIR}"
echo ""

for ((i=0; i<N_XYZ; i++)); do
    IDX=$(printf "%03d" $i)
    GEOM_DIR="${CURVE_DIR}/geom_${IDX}"
    XYZ="${XYZ_FILES[$i]}"
    RASORB="${RASORB_FILES[$i]}"

    mkdir -p "${GEOM_DIR}"

    # --- Per-geometry config ---
    # Patch xyz_file, rasorb_path, account; keep all other settings from config_Po2.yml
    sed \
        -e "s|xyz_file:.*|xyz_file: \"${XYZ}\"|" \
        -e "s|rasorb_path:.*|rasorb_path: \"${RASORB}\"|" \
        -e "s|account:.*|account: \"2025_060\"|" \
        "${WORKFLOW_DIR}/config_Po2.yml" > "${GEOM_DIR}/config_geom_${IDX}.yml"

    # --- Per-geometry submit script (generated fresh, not sed-patched) ---
    cat > "${GEOM_DIR}/submit_geom_${IDX}.sh" << SUBMITEOF
#!/bin/bash
#PBS -N Po2_EFFE_geom_${IDX}
#PBS -o submit_geom_${IDX}.out
#PBS -e submit_geom_${IDX}.err
#PBS -l walltime=02:30:00
#PBS -l mem=1gb
#PBS -l nodes=1:ppn=1
#PBS -m be -M joachim.scheerlinck@ugent.be
#PBS -A 2025_060

echo "=========================================="
echo "Po2 EFFE MS-CASPT2 geom_${IDX}"
echo "Job ID: \$PBS_JOBID"
echo "Start time: \$(date)"
echo "=========================================="

ulimit -s unlimited
cd ${GEOM_DIR}

module load Python/3.11.3-GCCcore-12.3.0

echo "Activating virtual environment..."
source ${WORKFLOW_DIR}/mscaspt2_venv/bin/activate

echo "Python: \$(python --version)"
echo "Parsl:  \$(python -c 'import parsl; print(parsl.__version__)' 2>/dev/null || echo 'NOT FOUND')"
echo ""

echo "Starting Po2 EFFE workflow for geom_${IDX}..."
python ${WORKFLOW_DIR}/run_mscaspt2_Po2.py config_geom_${IDX}.yml

EXIT_CODE=\$?

echo ""
echo "End time: \$(date)"
if [ \$EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "geom_${IDX} completed successfully!"
    echo "=========================================="
else
    echo "=========================================="
    echo "ERROR: geom_${IDX} failed (exit code \$EXIT_CODE)"
    echo "=========================================="
fi

deactivate
exit \$EXIT_CODE
SUBMITEOF
    chmod +x "${GEOM_DIR}/submit_geom_${IDX}.sh"

    echo "  geom_${IDX}: $(basename ${XYZ}) <-> $(basename ${RASORB})"
done

echo ""
echo "Setup complete. Next steps:"
echo "  1. Spot-check: grep 'xyz_file\|rasorb_path' ${CURVE_DIR}/geom_000/config_geom_000.yml"
echo "     (should show po2_000.xyz paired with system_24.RasOrb)"
echo "  2. Submit all:  bash ${WORKFLOW_DIR}/submit_all_po2.sh"
echo "  3. Or batch:    bash ${WORKFLOW_DIR}/submit_all_po2.sh 0 9"
