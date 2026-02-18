# MS-CASPT2 Root-Specific Workflow with Parsl

This workflow automates the process of running root-specific MS-CASPT2 calculations, extracting couplings, and combining results using the EFFE keyword.

## Files Overview

- `run_mscaspt2_workflow.py` - Main Parsl workflow script
- `config.yml` - Configuration file with all calculation parameters
- `submit_workflow.sh` - PBS submission script
- `sorbaldehyde.xyz` - Your molecule XYZ file

## Setup

### 1. Create Virtual Environment (One-Time Setup)

Since autoCAS module is incompatible with the Python/Parsl/NumPy modules, we use a virtual environment:

```bash
# Run the setup script
bash setup_venv.sh
```

This will:
- Create a virtual environment called `mscaspt2_venv`
- Install Parsl, NumPy, and PyYAML in the virtual environment
- The virtual environment will be automatically activated by `submit_workflow.sh`

**Important**: The autoCAS module will be loaded in the **worker jobs** (not the master job), so pymolcas will be available where it's needed.

### 2. Prepare Your Files

Make sure you have:
- RasOrb files from autoCAS for each spin state (singlet, triplet, quintet)
- XYZ file for your molecule
- RASSCF h5 files for the EFFE step

### 3. Edit Configuration

Edit `config.yml` to match your system:

```yaml
cluster:
  account: "YOUR_ACCOUNT"  # Your PBS account
  workdir_base: "/path/to/your/workdir"  # Your scratch directory

molecule:
  xyz_file: "your_molecule.xyz"

calculations:
  - name: "singlets"
    rasorb_path: "path/to/singlet.RasOrb"  # From autoCAS
    rasscf_h5: "singlet.rasscf.h5"
    n_roots: 3
    # ... other parameters
```

## Usage

### Submit the workflow:

```bash
qsub submit_workflow.sh
```

### What happens:

1. The master job starts on one node
2. Parsl reads your config and spawns worker jobs
3. For each spin state (singlet/triplet/quintet):
   - Creates input files for each root
   - Submits parallel MOLCAS jobs (one per root)
   - Extracts couplings from log files
   - Combines results with EFFE keyword
   - Creates final RASSI input
4. Output files are organized by spin state

## Directory Structure

After running, you'll have:

```
./
├── singlets/
│   ├── root1/
│   │   ├── root1.inp
│   │   ├── root1.log
│   │   └── root1.out
│   ├── root2/
│   ├── root3/
│   └── combined/
│       └── combined_mscaspt2.inp  # Final input with EFFE
├── triplets/
│   └── ...
└── quintets/
    └── ...
```

## Configuration Options

### Cluster Settings

- `account`: Your PBS account
- `walltime`: Time limit per worker job
- `max_blocks`: Maximum parallel jobs (adjust based on cluster load)
- `molcas_mem`: Memory per MOLCAS job (MB)
- `molcas_nprocs`: Processors per MOLCAS job

### Calculation Settings (per spin state)

- `name`: Identifier (e.g., "singlets", "triplets")
- `spin`: Spin multiplicity (1=singlet, 3=triplet, 5=quintet)
- `n_roots`: Number of roots to calculate
- `rasorb_path`: Path to RasOrb file from autoCAS
- `rasscf_h5`: RASSCF h5 file for EFFE step
- `inactive`: Number of inactive orbitals
- `ras2`: Number of RAS2 orbitals
- `nactel`: Number of active electrons
- `basis`: Basis set (default: ANO-S-VDZP)
- `imaginary`: Imaginary shift (default: 0.1)
- `output_dir`: Where to store results

## Monitoring

Check workflow progress:
```bash
# Check master job
qstat -u $USER

# View output
tail -f workflow.out

# Check individual root jobs (spawned by Parsl)
qstat -u $USER | grep parsl
```

## Troubleshooting

### "Module not found" errors
Make sure autoCAS and Python modules are available on your cluster.

### "Permission denied" on workdir
Check that `workdir_base` in config.yml points to a writable scratch directory.

### Jobs not spawning
- Verify your PBS account is correct
- Check cluster queue limits
- Try reducing `max_blocks` if cluster is busy

### Coupling extraction fails
- Check that MOLCAS jobs completed successfully
- Verify log files exist in `root*/root*.log`
- Ensure "Hamiltonian Effective Couplings" section is in log files

## Next Steps: SOCASSI

After completing this workflow, you'll have RASSI outputs for each spin state:
- `singlets/combined/JOB001` (JobIph file)
- `triplets/combined/JOB001`
- `quintets/combined/JOB001`

These can be combined in a final SOCASSI calculation to include spin-orbit coupling. Let me know if you need help with that step!

## Customization

To run only specific spin states, edit `config.yml` and comment out unwanted calculations:

```yaml
calculations:
  - name: "singlets"
    # ... singlet config
  
  # - name: "triplets"  # Commented out
  #   # ... triplet config
```

## Performance Tips

1. **Parallel efficiency**: Set `max_blocks` based on:
   - Number of roots × number of spin states
   - Cluster availability
   - Example: 3 singlets + 5 triplets + 2 quintets = 10 roots = set `max_blocks: 10`

2. **Memory**: Adjust `molcas_mem` based on system size:
   - Small systems: 2000-4000 MB
   - Medium systems: 4000-8000 MB
   - Large systems: 8000-16000 MB

3. **Walltime**: Root calculations typically take 5-30 minutes each. Set worker walltime accordingly.
