#!/bin/bash
#PBS -N Po2_EFFE_workflow
#PBS -o workflow_Po2.out
#PBS -e workflow_Po2.err
#PBS -l walltime=16:00:00
#PBS -l mem=1gb
#PBS -l nodes=1:ppn=1
#PBS -m be -M joachim.scheerlinck@ugent.be
#PBS -A starting_2025_097

# =============================================================================
# Master orchestrator for Po2 EFFE MS-CASPT2 workflow.
# Pipelined full-RASSCF mode (default):
#   Singlet RASSCF → submit 28 singlet CASPT2 jobs (SLURM, parallel)
#                  + Triplet RASSCF → submit 90 triplet CASPT2 jobs
#                  + Quintet RASSCF → submit 70 quintet CASPT2 jobs
#   → collect all + EFFE combine × 3 → final SO-RASSI
#
# walltime=48h: 3 RASSCF steps (~1h each) + wait for slowest CASPT2 block
#               (~90 roots × serial at most 2h = ~2h if all parallel) + overhead.
#
# BEFORE SUBMITTING:
#   1. Update rasorb_path entries in config_Po2.yml to ABSOLUTE paths on HPC
#   2. Verify xyz_file path in config_Po2.yml
#   3. Run: bash setup_venv.sh (if not already done)
# =============================================================================

echo "=========================================="
echo "Po2 EFFE MS-CASPT2 Workflow Starting"
echo "Job ID: $PBS_JOBID"
echo "Working directory: $PBS_O_WORKDIR"
echo "Start time: $(date)"
echo "=========================================="

ulimit -s unlimited
cd ${PBS_O_WORKDIR}

# Activate virtual environment (created by setup_venv.sh)
if [ ! -d "mscaspt2_venv" ]; then
    echo "ERROR: Virtual environment not found in ${PBS_O_WORKDIR}"
    echo "Please run: bash setup_venv.sh"
    exit 1
fi

module load Python/3.11.3-GCCcore-12.3.0

echo "Activating virtual environment..."
source mscaspt2_venv/bin/activate

# Verify
echo "Python: $(python --version)"
echo "Parsl:  $(python -c 'import parsl; print(parsl.__version__)' 2>/dev/null || echo 'NOT FOUND')"
echo ""

# Run the Po2 EFFE workflow
echo "Starting Po2 EFFE workflow..."
python run_mscaspt2_Po2.py config_Po2.yml

EXIT_CODE=$?

echo ""
echo "End time: $(date)"
if [ $EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "Workflow completed successfully!"
    echo "=========================================="
else
    echo "=========================================="
    echo "ERROR: Workflow failed (exit code $EXIT_CODE)"
    echo "=========================================="
fi

deactivate
exit $EXIT_CODE
