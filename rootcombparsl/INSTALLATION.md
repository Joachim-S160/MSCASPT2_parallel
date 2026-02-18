# Quick Installation Guide

## Problem: Module Incompatibility

The autoCAS module (`autoCAS/3.0.0-iomkl-2023a`) cannot be loaded simultaneously with:
- Parsl/2023.7.17-GCCcore-11.3.0
- numpy modules
- PyYAML modules

## Solution: Virtual Environment

We use a **virtual environment** for the master workflow job, and load autoCAS only in the worker jobs where pymolcas runs.

## Step-by-Step Installation

### 1. Download All Files

You need these files in your working directory:
- `run_mscaspt2_workflow.py` - Main workflow script
- `config.yml` - Configuration file
- `setup_venv.sh` - Virtual environment setup script
- `submit_workflow.sh` - PBS submission script
- `sorbaldehyde.xyz` - Your molecule file

### 2. Create Virtual Environment

```bash
# Make the setup script executable
chmod +x setup_venv.sh

# Run it (only needed once)
bash setup_venv.sh
```

Expected output:
```
Setting up virtual environment for MS-CASPT2 workflow...
Creating virtual environment...
Installing Parsl...
Installing numpy...
==========================================
Virtual environment setup complete!
==========================================
```

### 3. Verify Installation

```bash
# Activate the virtual environment manually to test
source mscaspt2_venv/bin/activate

# Check installations
python -c "import parsl; print('Parsl version:', parsl.__version__)"
python -c "import numpy; print('NumPy version:', numpy.__version__)"
python -c "import yaml; print('PyYAML installed')"

# Deactivate
deactivate
```

### 4. Edit Configuration

Edit `config.yml`:

```yaml
cluster:
  account: "2025_026"  # Your account
  module: "autoCAS/3.0.0-iomkl-2023a"  # Keep this exact name
  workdir_base: "/path/to/your/scratch"

molecule:
  xyz_file: "sorbaldehyde.xyz"

calculations:
  - name: "singlets"
    rasorb_path: "/full/path/to/sorbaldehyde-singlet.RasOrb"
    rasscf_h5: "singlet.rasscf.h5"
    # ... adjust other parameters
```

**Important**: Use **full paths** for `rasorb_path` since worker jobs run in different directories!

### 5. Submit Workflow

```bash
qsub submit_workflow.sh
```

### 6. Monitor Progress

```bash
# Check master job
qstat -u $USER

# Watch output
tail -f workflow.out

# Check for worker jobs
watch -n 5 'qstat -u $USER'
```

## How It Works

### Master Job (submit_workflow.sh)
- Runs on 1 node
- Uses virtual environment (Python + Parsl + NumPy + PyYAML)
- **Does NOT load autoCAS** (to avoid conflicts)
- Coordinates the workflow

### Worker Jobs (spawned by Parsl)
- Automatically submitted by Parsl for each root calculation
- **Load autoCAS module** (for pymolcas)
- Run individual MOLCAS calculations
- No module conflicts because each worker is independent

## Troubleshooting

### "Virtual environment not found"
```bash
# Run setup again
bash setup_venv.sh
```

### "pymolcas not found" in worker jobs
Check that `config.yml` has the correct autoCAS module:
```yaml
cluster:
  module: "autoCAS/3.0.0-iomkl-2023a"
```

### Worker jobs fail immediately
Check `workflow.err` and look for PBS errors. Common issues:
- Wrong account name
- Insufficient walltime
- Queue restrictions

### Python import errors in master job
Verify virtual environment:
```bash
source mscaspt2_venv/bin/activate
pip list | grep -E "parsl|numpy|PyYAML"
deactivate
```

If packages are missing:
```bash
source mscaspt2_venv/bin/activate
pip install parsl numpy pyyaml
deactivate
```

## File Checklist

Before submitting:

- [ ] Virtual environment created (`mscaspt2_venv/` directory exists)
- [ ] `config.yml` edited with your paths and settings
- [ ] XYZ file exists and path is correct in config
- [ ] RasOrb files exist (from autoCAS) - use full paths!
- [ ] RASSCF h5 files exist or will be generated
- [ ] PBS account is correct
- [ ] Scratch directory (`workdir_base`) is writable

## Expected Runtime

For typical systems:
- Setup: 2-3 minutes (one time only)
- Master job: Runs for duration of workflow
- Worker jobs: 5-30 minutes per root
- Total: Depends on number of roots and parallel capacity

Example: 3 singlets + 5 triplets + 2 quintets = 10 roots
- Sequential: ~3 hours
- Parallel (max_blocks=10): ~30 minutes

## Clean Up

After successful runs, you can clean up working directories:
```bash
# Remove MOLCAS workdirs (they can be large!)
rm -rf /dodrio/scratch/projects/starting_2025_097/molcas_workdirs/parsl_*
```

Keep your virtual environment for future runs!
