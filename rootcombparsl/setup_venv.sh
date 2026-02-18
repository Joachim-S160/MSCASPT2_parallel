#!/bin/bash
# Setup script to create virtual environment with required packages
# Run this ONCE before submitting your first workflow

echo "Setting up virtual environment for MS-CASPT2 workflow..."

# Load compatible modules for setup
module load Python/3.11.3-GCCcore-12.3.0
module load PyYAML/6.0-GCCcore-12.3.0
#module load SciPy-bundle/2025.06-gfbf-2025a
# Create virtual environment
echo "Creating virtual environment..."
python -m venv mscaspt2_venv

# Activate it
source mscaspt2_venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required packages
echo "Installing Parsl..."
pip install parsl

echo "Installing numpy..."
pip install numpy

echo "Installing PyYAML..."
pip install pyyaml

echo ""
echo "=========================================="
echo "Virtual environment setup complete!"
echo "=========================================="
echo ""
echo "The virtual environment is located at: ./mscaspt2_venv"
echo ""
echo "To use it in your workflow, it will be automatically activated"
echo "by the submit_workflow.sh script."
echo ""
echo "You can now submit your workflow with: qsub submit_workflow.sh"
