#!/usr/bin/env python3
"""
Parsl workflow for root-specific MS-CASPT2 calculations
Automates splitting roots, running individual CASPT2, extracting couplings, and recombining
"""
import parsl
from parsl import python_app, bash_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider
from parsl.launchers import MpiRunLauncher
import os
import re
import yaml
from pathlib import Path
from typing import List, Tuple, Dict, Any
import numpy as np
import argparse


def load_config(config_file: str) -> Dict[str, Any]:
    """Load configuration from YAML file"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def setup_parsl(cluster_config: Dict[str, Any]) -> None:
    """Configure Parsl based on cluster settings"""
    config = Config(
        executors=[
            HighThroughputExecutor(
                label="molcas_executor",
                cores_per_worker=1,  # Each worker gets 1 core
                provider=SlurmProvider(
                    partition=cluster_config.get('partition', 'cpu_milan_rhel9'),
                    account=cluster_config['account'],
                    nodes_per_block=1,
                    init_blocks=cluster_config.get('init_blocks', 0),
                    max_blocks=cluster_config.get('max_blocks', 20),
                    walltime=cluster_config.get('walltime', '00:30:00'),
                    scheduler_options=cluster_config.get(
                        'scheduler_options',
                        "#SBATCH --clusters=dodrio\n"
                    ),
                    worker_init=f"""
# Load autoCAS module (required for pymolcas)
module purge
module load {cluster_config.get('module', 'autoCAS/3.0.0-iomkl-2023a')}

# Set MOLCAS environment
export MOLCAS_PRINT={cluster_config.get('molcas_print', 'verbose')}
export MOLCAS_MEM={cluster_config.get('molcas_mem', 4000)}
export MOLCAS_NPROCS={cluster_config.get('molcas_nprocs', 1)}
export OMP_NUM_THREADS={cluster_config.get('omp_threads', 1)}
ulimit -s unlimited

# Verify pymolcas is available
which pymolcas || echo "WARNING: pymolcas not found in PATH"
                    """,
                    launcher=MpiRunLauncher(),
                ),
            )
        ],
        strategy='htex_auto_scale',
    )
    parsl.load(config)


@python_app
def create_root_input(
    root_idx: int,
    calc_params: Dict[str, Any],
    output_dir: str
) -> str:
    """Create input file for a specific root MS-CASPT2 calculation"""
    from pathlib import Path
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    input_file = f"{output_dir}/root{root_idx}.inp"
    xyz_file = f"{output_dir}/molecule.xyz"
    
    # Write XYZ file
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])
    import os
    # Use absolute paths
    output_dir = os.path.abspath(output_dir)
    # Use absolute path for RasOrb
    rasorb_path = os.path.abspath(calc_params['rasorb_path'])
    
    # Create input file
    input_content = f"""&GATEWAY
Title
MS-CASPT2 Root {root_idx} of {calc_params['n_roots']} (Spin={calc_params['spin']})
COORD = {xyz_file}
GROUP = nosymm
BASIS = {calc_params['basis']}
RICD
&SEWARD
DoAnalytical
&RASSCF
File = {rasorb_path}
CIONLY
SPIN={calc_params['spin']}
CIROOT = {calc_params['n_roots']} {calc_params['n_roots']} 1
Inactive = {calc_params['inactive']}
RAS2 = {calc_params['ras2']}
NACTEL = {calc_params['nactel']}
ORBAppear = COMPACT
&CASPT2
Multistate = all
Imaginary = {calc_params['imaginary']}
only={root_idx}
"""
    
    with open(input_file, 'w') as f:
        f.write(input_content)
    
    return input_file


@bash_app
def run_molcas(input_file: str, workdir_base: str, nprocs: int, inputs=[]) -> str:
    """Run MOLCAS calculation for a single root"""
    import os
    from pathlib import Path
    
    input_path = Path(input_file)
    output_dir = input_path.parent
    job_name = input_path.stem
    
    # Create unique workdir
    import time
    unique_id = f"parsl_{os.getpid()}_{int(time.time())}"
    workdir = f"{workdir_base}/{unique_id}_{job_name}"
    
    return f"""
    set -e
    export MOLCAS_WORKDIR={workdir}
    mkdir -p $MOLCAS_WORKDIR
    cd {output_dir}
    pymolcas -np {nprocs} {input_file} -f
    echo {output_dir}/{job_name}.log
    """


@python_app
def extract_couplings(log_file: str, inputs=[]) -> Tuple[int, List[float]]:
    """Extract effective couplings from MOLCAS log file"""
    import re
    from pathlib import Path
    
    log_path = Path(log_file)
    root_match = re.search(r'root(\d+)', log_path.stem)
    if not root_match:
        raise ValueError(f"Cannot extract root index from {log_path.stem}")
    root_idx = int(root_match.group(1))
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Find the Hamiltonian Effective Couplings section
    pattern = r'Hamiltonian Effective Couplings.*?\n((?:.*?<.*?\n)*)'
    match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
    
    if not match:
        raise ValueError(f"Could not find effective couplings in {log_file}")
    
    # Extract coupling values
    coupling_lines = match.group(1)
    couplings = []
    for line in coupling_lines.split('\n'):
        if '<' in line:
            # Extract the last column (coupling value)
            value = float(line.split()[-1])
            couplings.append(value)
    
    return (root_idx, couplings)


@python_app
def combine_roots(
    coupling_data: List[Tuple[int, List[float]]],
    calc_params: Dict[str, Any],
    output_dir: str,
    inputs=[]
) -> str:
    """Create combined MS-CASPT2 input with EFFE keyword"""
    from pathlib import Path
    import numpy as np
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Sort by root index
    coupling_data = sorted(coupling_data, key=lambda x: x[0])
    
    n_roots = calc_params['n_roots']
    
    # Build effective Hamiltonian matrix
    H_eff = np.zeros((n_roots, n_roots))
    for root_idx, couplings in coupling_data:
        if len(couplings) != n_roots:
            raise ValueError(f"Root {root_idx} has {len(couplings)} couplings, expected {n_roots}")
        H_eff[root_idx - 1, :] = couplings
    
    # Create EFFE section
    effe_lines = [f"{n_roots}"]
    for i in range(n_roots):
        row_str = " ".join([f"{H_eff[i, j]:.14E}" for j in range(n_roots)])
        effe_lines.append(row_str)
    effe_block = "\n".join(effe_lines)
    
    # Create combined input
    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    input_file = f"{output_dir}/combined_mscaspt2.inp"
    input_content = f"""&GATEWAY
Title
Combined MS-CASPT2 with EFFE (Spin={calc_params['spin']})
COORD = {xyz_file}
GROUP = nosymm
BASIS = {calc_params['basis']}
RICD
&SEWARD
DoAnalytical
&CASPT2
File={calc_params['rasscf_h5']}
Multistate = all
Imaginary = {calc_params['imaginary']}
EFFE
{effe_block}
>> COPY $Project.caspt2.h5 JOB001
&RASSI
NR OF JOBIPHS = 1 all
TRDI
"""
    
    with open(input_file, 'w') as f:
        f.write(input_content)
    
    return input_file


def process_spin_state(
    spin_config: Dict[str, Any],
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int
) -> str:
    """Process a single spin state (singlet, triplet, or quintet)"""
    print(f"Processing spin state: {spin_config['name']} (spin={spin_config['spin']})")
    
    # Prepare calculation parameters
    calc_params = {
        'rasorb_path': spin_config['rasorb_path'],
        'xyz_content': xyz_content,
        'n_roots': spin_config['n_roots'],
        'spin': spin_config['spin'],
        'inactive': spin_config['inactive'],
        'ras2': spin_config['ras2'],
        'nactel': spin_config['nactel'],
        'basis': spin_config.get('basis', 'ANO-S-VDZP'),
        'imaginary': spin_config.get('imaginary', 0.1),
        'rasscf_h5': spin_config['rasscf_h5']
    }
    
    base_output_dir = spin_config['output_dir']
    
    # Create root-specific input files
    print(f"  Creating {calc_params['n_roots']} root input files...")
    input_futures = []
    for root_idx in range(1, calc_params['n_roots'] + 1):
        root_dir = f"{base_output_dir}/root{root_idx}"
        future = create_root_input(
            root_idx=root_idx,
            calc_params=calc_params,
            output_dir=root_dir
        )
        input_futures.append(future)
    
    # Wait for input files to be created
    input_files = [f.result() for f in input_futures]
    print(f"  Input files created: {len(input_files)}")
    
    # Run MOLCAS calculations in parallel
    print(f"  Submitting {len(input_files)} MOLCAS jobs...")
    molcas_futures = []
    for input_file in input_files:
        future = run_molcas(input_file, workdir_base, molcas_nprocs, inputs=[])
        molcas_futures.append(future)
    
    # Extract log file paths
    print(f"  Waiting for MOLCAS jobs to complete...")
    print(f"\n Futures: {[f for f in molcas_futures]} \n")
    log_files = [f.result().strip() for f in molcas_futures]
    print(f"  All jobs completed. Log files: {len(log_files)}")
    
    # Extract couplings from all roots
    print(f"  Extracting couplings...")
    coupling_futures = []
    for log_file in log_files:
        future = extract_couplings(log_file, inputs=[])
        coupling_futures.append(future)
    
    # Wait for all coupling extractions
    coupling_data = [f.result() for f in coupling_futures]
    print(f"  Couplings extracted: {len(coupling_data)}")
    
    # Combine results
    print(f"  Combining roots with EFFE...")
    combine_dir = f"{base_output_dir}/combined"
    combined_future = combine_roots(
        coupling_data=coupling_data,
        calc_params=calc_params,
        output_dir=combine_dir,
        inputs=[]
    )
    
    combined_input = combined_future.result()
    print(f"  Combined input created: {combined_input}")
    
    return combined_input


def main():
    parser = argparse.ArgumentParser(description='Run root-specific MS-CASPT2 workflow')
    parser.add_argument('config', help='Path to YAML configuration file')
    args = parser.parse_args()
    
    # Load configuration
    print(f"Loading configuration from {args.config}...")
    config = load_config(args.config)
    
    # Setup Parsl
    print("Setting up Parsl...")
    setup_parsl(config['cluster'])
    
    # Load XYZ file
    print(f"Loading XYZ file: {config['molecule']['xyz_file']}")
    with open(config['molecule']['xyz_file'], 'r') as f:
        xyz_content = f.read()
    
    # Process each spin state
    results = {}
    for spin_state in config['calculations']:
        combined_input = process_spin_state(
            spin_config=spin_state,
            xyz_content=xyz_content,
            workdir_base=config['cluster']['workdir_base'],
            molcas_nprocs=config['cluster'].get('molcas_nprocs', 1)
        )
        results[spin_state['name']] = combined_input
    
    # Summary
    print("\n" + "="*60)
    print("WORKFLOW COMPLETED SUCCESSFULLY")
    print("="*60)
    for name, path in results.items():
        print(f"{name}: {path}")
    print("="*60)
    
    # Cleanup Parsl
    parsl.clear()


if __name__ == "__main__":
    main()
