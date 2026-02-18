#!/bin/bash
#PBS -N MSCASPT2_workflow
#PBS -o workflow.out
#PBS -e workflow.err
#PBS -l walltime=02:00:00
#PBS -l mem=8gb
#PBS -m be -M joachim.scheerlinck@ugent.be
#PBS -l nodes=1:ppn=1
#PBS -A 2025_026

# This is the master job that coordinates the workflow
# Individual MOLCAS jobs will be spawned by Parsl

echo "=========================================="
echo "MS-CASPT2 Workflow Starting"
echo "Job ID: $PBS_JOBID"
echo "Working directory: $PBS_O_WORKDIR"
echo "=========================================="

# Set unlimited stack size
ulimit -s unlimited

# Change to working directory
cd ${PBS_O_WORKDIR}

module load Python/3.11.3-GCCcore-12.3.0
module load PyYAML/6.0-GCCcore-12.3.0

echo "Activating virtual environment..."
source /dodrio/scratch/projects/starting_2025_097/speedup_mscaspt2/rootcombparsl/mscaspt2_venv/bin/activate

# Verify installations
echo "Python version: $(python --version)"
echo "Parsl installed: $(python -c 'import parsl; print(parsl.__version__)' 2>/dev/null || echo 'NOT FOUND')"
echo "NumPy installed: $(python -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo 'NOT FOUND')"
echo "PyYAML installed: $(python -c 'import yaml; print(yaml.__version__)' 2>/dev/null || echo 'NOT FOUND')"
echo ""

# Run the workflow
echo "Starting workflow execution..."
python simple_parallel_molcas.py

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Workflow completed successfully!"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "ERROR: Workflow failed!"
    echo "=========================================="
    exit 1
fi

# Deactivate virtual environment
deactivate
