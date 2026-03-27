#!/usr/bin/env python3
"""
Po2 MS-CASPT2 EFFE workflow — 6-block C2 symmetry SO-CASSI.

Adapts run_mscaspt2_workflow.py for Po2's 6-spin/symmetry-block structure:
  - 6 blocks: singlet irr1/2, triplet irr1/2, quintet irr1/2
  - C2 symmetry (Group=XY), two-irrep Inactive/Ras2 format
  - Binary JOBMIX output (HPC OpenMolcas without HDF5 mode)
  - Final SO-RASSI combining JOB001-JOB006 WITHOUT Ejob keyword

Two modes (--full-rasscf flag):
  Default (CIONLY): uses pre-converged autoCAS orbitals directly for CASPT2.
  Full RASSCF:      runs state-specific RASSCF per spin group first, then passes
                    optimized orbitals to CASPT2. Orbital handoff between spin
                    states: singlet → triplet → quintet (each uses previous as guess).

Workflow per block:
  n_roots individual CASPT2(only=N) jobs [parallel via Parsl SLURM]
  → extract "Hamiltonian Effective Couplings" from each log
  → combined CASPT2+EFFE job → >>COPY $Project.JobMix JOB00X
Final step:
  SO-RASSI over JOB001-JOB006 (Spin Orbit, Omega, no Ejob)

Config: config_Po2.yml (6-block structure with inactive_c2, ras2_c2, job_number per block)

Fixes vs. original:
  - MpiRunLauncher → SimpleLauncher (single-node jobs, no MPI)
  - Log path: pymolcas writes log to input dir, not MOLCAS_WORKDIR
  - Added full RASSCF + orbital handoff option (--full-rasscf)
"""

import parsl
from parsl import python_app, bash_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider
from parsl.launchers import SimpleLauncher
import os
import yaml
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import numpy as np
import argparse


def load_config(config_file: str) -> Dict[str, Any]:
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def setup_parsl(cluster_config: Dict[str, Any]) -> None:
    """Configure Parsl with SLURM provider for HPC."""
    config = Config(
        executors=[
            HighThroughputExecutor(
                label="molcas_executor",
                cores_per_worker=1,
                provider=SlurmProvider(
                    partition=cluster_config.get('partition', 'cpu_milan_rhel9'),
                    account=cluster_config['account'],
                    nodes_per_block=1,
                    init_blocks=cluster_config.get('init_blocks', 0),
                    max_blocks=cluster_config.get('max_blocks', 220),
                    walltime=cluster_config.get('walltime', '02:00:00'),
                    scheduler_options=cluster_config.get(
                        'scheduler_options', "#SBATCH --clusters=dodrio\n"
                    ),
                    worker_init=f"""
module purge
module load {cluster_config.get('module', 'autoCAS/3.0.0-iomkl-2023a')}
export MOLCAS_PRINT={cluster_config.get('molcas_print', 'verbose')}
export MOLCAS_MEM={cluster_config.get('molcas_mem', 20000)}
export MOLCAS_NPROCS={cluster_config.get('molcas_nprocs', 1)}
export OMP_NUM_THREADS={cluster_config.get('omp_threads', 1)}
ulimit -s unlimited
which pymolcas || echo "WARNING: pymolcas not found in PATH"
""",
                    launcher=SimpleLauncher(),  # Fixed: was MpiRunLauncher()
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
    output_dir: str,
) -> str:
    """Create RASSCF(CIONLY)+CASPT2(only=N) input for one root of one block."""
    from pathlib import Path
    import os

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    rasorb_path = os.path.abspath(calc_params['rasorb_path'])

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    input_content = f"""&GATEWAY
Title
Po2 Root {root_idx}/{calc_params['n_roots']} spin={calc_params['spin']} sym={calc_params['symmetry']}
COORD = {xyz_file}
GROUP = XY
BASIS = {calc_params['basis']}
RICD
&SEWARD
Cholesky
&RASSCF
File = {rasorb_path}
CIONLY
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {calc_params['n_roots']} {calc_params['n_roots']} 1
Inactive = {calc_params['inactive_c2']}
Ras2 = {calc_params['ras2_c2']}
NACTEL = {calc_params['nactel']}
ORBAppear = COMPACT
&CASPT2
MAXITER = 300
Frozen = {calc_params['inactive_c2']}
Multistate = all
Imaginary = {calc_params['imaginary']}
only = {root_idx}
"""

    input_file = f"{output_dir}/root{root_idx}.inp"
    with open(input_file, 'w') as f:
        f.write(input_content)
    return input_file


@python_app
def create_rasscf_input(
    calc_params: Dict[str, Any],
    start_rasorb: str,
    output_dir: str,
    inputs=[],
) -> str:
    """Create full RASSCF input (no CIONLY) for orbital optimization of one spin/irrep block."""
    from pathlib import Path
    import os

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    rasorb_path = os.path.abspath(start_rasorb)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    input_content = f"""&GATEWAY
Title
Po2 RASSCF spin={calc_params['spin']} sym={calc_params['symmetry']} (full opt)
COORD = {xyz_file}
GROUP = XY
BASIS = {calc_params['basis']}
RICD
&SEWARD
Cholesky
&RASSCF
File = {rasorb_path}
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {calc_params['n_roots']} {calc_params['n_roots']} 1
Inactive = {calc_params['inactive_c2']}
Ras2 = {calc_params['ras2_c2']}
NACTEL = {calc_params['nactel']}
Levshft = 0.5
Iteration = 200 50
ORBAppear = COMPACT
"""

    input_file = f"{output_dir}/rasscf_full.inp"
    with open(input_file, 'w') as f:
        f.write(input_content)
    return input_file


@bash_app
def run_molcas(input_file: str, workdir_base: str, nprocs: int, inputs=[]) -> str:
    """Run MOLCAS for one input file. Returns path to log file."""
    import os
    import time
    from pathlib import Path

    input_path = Path(input_file)
    output_dir = str(input_path.parent)
    job_name = input_path.stem
    unique_id = f"parsl_{os.getpid()}_{int(time.time())}"
    workdir = f"{workdir_base}/{unique_id}_{job_name}"

    # Fixed: log is written by pymolcas to the input dir, not MOLCAS_WORKDIR
    return f"""
set -e
export MOLCAS_WORKDIR={workdir}
mkdir -p $MOLCAS_WORKDIR
cd {output_dir}
pymolcas -np {nprocs} {input_file} -f
echo {output_dir}/{job_name}.log
"""


@bash_app
def run_rasscf_and_copy_orb(input_file: str, workdir_base: str, nprocs: int,
                             out_rasorb: str, inputs=[]) -> str:
    """Run full RASSCF and copy the output .RasOrb to a stable path.

    Returns the path to the copied .RasOrb file (out_rasorb).
    """
    import os
    import time
    from pathlib import Path

    input_path = Path(input_file)
    output_dir = str(input_path.parent)
    job_name = input_path.stem
    unique_id = f"parsl_{os.getpid()}_{int(time.time())}"
    workdir = f"{workdir_base}/{unique_id}_{job_name}"

    return f"""
set -e
export MOLCAS_WORKDIR={workdir}
mkdir -p $MOLCAS_WORKDIR
cd {output_dir}
pymolcas -np {nprocs} {input_file} -f
# pymolcas writes {job_name}.RasOrb in the current directory
if [ ! -f {output_dir}/{job_name}.RasOrb ]; then
    echo "ERROR: {output_dir}/{job_name}.RasOrb not found after RASSCF"
    ls -la {output_dir}/
    exit 1
fi
cp {output_dir}/{job_name}.RasOrb {out_rasorb}
echo {out_rasorb}
"""


@python_app
def extract_couplings(log_file: str, inputs=[]) -> Tuple[int, List[float]]:
    """Parse 'Hamiltonian Effective Couplings' from a CASPT2(only=N) log file."""
    import re
    from pathlib import Path

    log_path = Path(log_file)
    root_match = re.search(r'root(\d+)', log_path.stem)
    if not root_match:
        raise ValueError(f"Cannot extract root index from log stem: {log_path.stem}")
    root_idx = int(root_match.group(1))

    with open(log_file, 'r') as f:
        content = f.read()

    # Section header, then lines of the form:  < N | <value>
    pattern = r'Hamiltonian Effective Couplings.*?\n((?:.*?<.*?\n)*)'
    match = re.search(pattern, content, re.DOTALL | re.IGNORECASE)
    if not match:
        raise ValueError(f"Could not find 'Hamiltonian Effective Couplings' in {log_file}")

    couplings = []
    for line in match.group(1).split('\n'):
        if '<' in line:
            couplings.append(float(line.split()[-1]))

    return (root_idx, couplings)


@python_app
def create_combined_input(
    coupling_data: List[Tuple[int, List[float]]],
    calc_params: Dict[str, Any],
    output_dir: str,
    inputs=[],
) -> str:
    """Create RASSCF(CIONLY)+CASPT2(EFFE) combined input for one block.

    Assembles the full H_eff matrix from individual root couplings and embeds
    it as an EFFE block. Writes binary JOBMIX to $CurrDir/JOBxxx.
    """
    from pathlib import Path
    import numpy as np
    import os

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    rasorb_path = os.path.abspath(calc_params['rasorb_path'])

    coupling_data = sorted(coupling_data, key=lambda x: x[0])
    n_roots = calc_params['n_roots']

    H_eff = np.zeros((n_roots, n_roots))
    for root_idx, couplings in coupling_data:
        if len(couplings) != n_roots:
            raise ValueError(
                f"Root {root_idx}: got {len(couplings)} couplings, expected {n_roots}"
            )
        H_eff[root_idx - 1, :] = couplings

    effe_lines = [str(n_roots)]
    for i in range(n_roots):
        effe_lines.append(
            " ".join([f"{H_eff[i, j]:.14E}" for j in range(n_roots)])
        )
    effe_block = "\n".join(effe_lines)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    job_num = calc_params['job_number']
    job_label = f"JOB{job_num:03d}"

    input_content = f"""&GATEWAY
Title
Po2 Combined EFFE spin={calc_params['spin']} sym={calc_params['symmetry']} -> {job_label}
COORD = {xyz_file}
GROUP = XY
BASIS = {calc_params['basis']}
RICD
&SEWARD
Cholesky
&RASSCF
File = {rasorb_path}
CIONLY
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {n_roots} {n_roots} 1
Inactive = {calc_params['inactive_c2']}
Ras2 = {calc_params['ras2_c2']}
NACTEL = {calc_params['nactel']}
ORBAppear = COMPACT
&CASPT2
MAXITER = 300
Frozen = {calc_params['inactive_c2']}
Multistate = all
Imaginary = {calc_params['imaginary']}
EFFE
{effe_block}
>>COPY $Project.JobMix $CurrDir/{job_label}
"""

    input_file = f"{output_dir}/combined_effe.inp"
    with open(input_file, 'w') as f:
        f.write(input_content)
    return input_file


@python_app
def create_final_rassi_input(
    block_results: List[Dict[str, Any]],
    output_dir: str,
    xyz_content: str,
    basis: str,
    inputs=[],
) -> str:
    """Create the final SO-RASSI input combining all 6 JOBxxx files (no Ejob)."""
    from pathlib import Path
    import os

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)

    block_results = sorted(block_results, key=lambda x: x['job_number'])
    n_blocks = len(block_results)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(xyz_content)

    copy_lines = []
    for blk in block_results:
        job_label = f"JOB{blk['job_number']:03d}"
        src = f"{blk['combined_dir']}/{job_label}"
        copy_lines.append(f">>COPY {src} {job_label}")
    copy_block = "\n".join(copy_lines)

    input_content = f"""&GATEWAY
Title
Po2 Final SO-RASSI {n_blocks} blocks (no Ejob)
COORD = {xyz_file}
GROUP = XY
BASIS = {basis}
RICD
&SEWARD
Cholesky
{copy_block}
&RASSI
Nr of JobIphs = {n_blocks} all
Spin Orbit
Omega
End of Input
"""

    input_file = f"{output_dir}/final_rassi.inp"
    with open(input_file, 'w') as f:
        f.write(input_content)
    return input_file


def process_block(
    spin_config: Dict[str, Any],
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int,
    rasorb_override: Optional[str] = None,
) -> Dict[str, Any]:
    """Process one spin/symmetry block: parallel root jobs → extract → combine.

    rasorb_override: if given, use this .RasOrb instead of spin_config['rasorb_path'].
    Returns dict with block metadata needed by the final RASSI step.
    """
    name = spin_config['name']
    n_roots = spin_config['n_roots']
    job_num = spin_config['job_number']
    rasorb_path = rasorb_override if rasorb_override else spin_config['rasorb_path']
    print(
        f"\n--- Block: {name} | "
        f"spin={spin_config['spin']}, sym={spin_config['symmetry']}, "
        f"n_roots={n_roots}, -> JOB{job_num:03d} ---"
    )
    if rasorb_override:
        print(f"  Using orbital override: {rasorb_override}")

    calc_params = {
        'rasorb_path': rasorb_path,
        'xyz_content': xyz_content,
        'n_roots': n_roots,
        'spin': spin_config['spin'],
        'symmetry': spin_config['symmetry'],
        'inactive_c2': spin_config['inactive_c2'],
        'ras2_c2': spin_config['ras2_c2'],
        'nactel': spin_config['nactel'],
        'basis': spin_config.get('basis', 'ANO-RCC-VQZP'),
        'imaginary': spin_config.get('imaginary', 0.25),
        'job_number': job_num,
    }

    base_dir = spin_config['output_dir']

    print(f"  Creating {n_roots} root input files...")
    input_futures = [
        create_root_input(root_idx, calc_params, f"{base_dir}/root{root_idx}")
        for root_idx in range(1, n_roots + 1)
    ]
    input_files = [f.result() for f in input_futures]
    print(f"  Input files ready ({len(input_files)} files).")

    print(f"  Submitting {n_roots} CASPT2(only=N) jobs to SLURM...")
    mol_futures = [
        run_molcas(inp, workdir_base, molcas_nprocs, inputs=[])
        for inp in input_files
    ]
    log_files = [f.result().strip() for f in mol_futures]
    print(f"  All {len(log_files)} root jobs completed.")

    print(f"  Extracting effective Hamiltonian couplings...")
    coupling_futures = [
        extract_couplings(log, inputs=[]) for log in log_files
    ]
    coupling_data = [f.result() for f in coupling_futures]
    print(f"  Couplings extracted ({len(coupling_data)} rows).")

    combine_dir = os.path.abspath(f"{base_dir}/combined")
    combined_inp_future = create_combined_input(
        coupling_data=coupling_data,
        calc_params=calc_params,
        output_dir=combine_dir,
        inputs=[],
    )
    combined_inp = combined_inp_future.result()
    print(f"  Combined EFFE input: {combined_inp}")

    print(f"  Running combined CASPT2+EFFE job...")
    combined_log_future = run_molcas(
        combined_inp, workdir_base, molcas_nprocs, inputs=[]
    )
    combined_log = combined_log_future.result().strip()
    print(f"  Combined job done: {combined_log}")

    return {
        'name': name,
        'job_number': job_num,
        'n_roots': n_roots,
        'combined_dir': combine_dir,
        'combined_log': combined_log,
    }


def process_spin_group_with_rasscf(
    blocks: List[Dict[str, Any]],
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int,
    start_rasorb: str,
) -> Tuple[List[Dict[str, Any]], str]:
    """Process all blocks in one spin group with full RASSCF + orbital handoff.

    Runs RASSCF for each irrep block sequentially (irr2 uses irr1 orbitals as
    starting guess), then runs all CASPT2 root jobs for irr1 and irr2 in parallel.

    Returns: (list of block results, path to last RASSCF .RasOrb for next spin group)
    """
    spin = blocks[0]['spin']
    print(f"\n=== Spin group S={spin}: full RASSCF + parallel CASPT2 ===")

    # Step 1: Run RASSCF for each irrep sequentially, passing orbitals forward
    rasorb_for_block = {}
    current_rasorb = start_rasorb
    for block in blocks:
        name = block['name']
        rasscf_dir = os.path.abspath(f"{block['output_dir']}/rasscf")

        calc_params = {
            'xyz_content': xyz_content,
            'spin': block['spin'],
            'symmetry': block['symmetry'],
            'n_roots': block['n_roots'],
            'inactive_c2': block['inactive_c2'],
            'ras2_c2': block['ras2_c2'],
            'nactel': block['nactel'],
            'basis': block.get('basis', 'ANO-RCC-VQZP'),
            'imaginary': block.get('imaginary', 0.25),
        }

        print(f"  RASSCF: {name} (spin={block['spin']}, sym={block['symmetry']}) "
              f"starting from {current_rasorb}")
        inp_future = create_rasscf_input(
            calc_params=calc_params,
            start_rasorb=current_rasorb,
            output_dir=rasscf_dir,
            inputs=[],
        )
        rasscf_inp = inp_future.result()

        out_rasorb = os.path.abspath(f"{rasscf_dir}/{name}.RasOrb")
        rasorb_log_future = run_rasscf_and_copy_orb(
            input_file=rasscf_inp,
            workdir_base=workdir_base,
            nprocs=molcas_nprocs,
            out_rasorb=out_rasorb,
            inputs=[],
        )
        rasorb_path = rasorb_log_future.result().strip()
        print(f"  RASSCF done: {rasorb_path}")

        rasorb_for_block[name] = rasorb_path
        current_rasorb = rasorb_path  # irr2 starts from irr1 output

    last_rasorb = current_rasorb  # return to caller for next spin group

    # Step 2: Run CASPT2 for all blocks in parallel (irr1 and irr2 simultaneously)
    print(f"  Submitting CASPT2 root jobs for all {len(blocks)} irrep blocks in parallel...")
    block_futures = {}
    for block in blocks:
        name = block['name']
        block_futures[name] = _submit_caspt2_for_block(
            spin_config=block,
            xyz_content=xyz_content,
            workdir_base=workdir_base,
            molcas_nprocs=molcas_nprocs,
            rasorb_path=rasorb_for_block[name],
        )

    # Step 3: Collect results
    block_results = []
    for block in blocks:
        name = block['name']
        result = block_futures[name]
        block_results.append(result)
        print(f"  Block {name} complete: JOB{result['job_number']:03d}")

    return block_results, last_rasorb


def _submit_caspt2_for_block(
    spin_config: Dict[str, Any],
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int,
    rasorb_path: str,
) -> Dict[str, Any]:
    """Submit all CASPT2 root jobs for one block and run EFFE combine. Returns block result."""
    name = spin_config['name']
    n_roots = spin_config['n_roots']
    job_num = spin_config['job_number']
    base_dir = spin_config['output_dir']

    calc_params = {
        'rasorb_path': rasorb_path,
        'xyz_content': xyz_content,
        'n_roots': n_roots,
        'spin': spin_config['spin'],
        'symmetry': spin_config['symmetry'],
        'inactive_c2': spin_config['inactive_c2'],
        'ras2_c2': spin_config['ras2_c2'],
        'nactel': spin_config['nactel'],
        'basis': spin_config.get('basis', 'ANO-RCC-VQZP'),
        'imaginary': spin_config.get('imaginary', 0.25),
        'job_number': job_num,
    }

    input_futures = [
        create_root_input(root_idx, calc_params, f"{base_dir}/root{root_idx}")
        for root_idx in range(1, n_roots + 1)
    ]
    input_files = [f.result() for f in input_futures]

    mol_futures = [
        run_molcas(inp, workdir_base, molcas_nprocs, inputs=[])
        for inp in input_files
    ]
    log_files = [f.result().strip() for f in mol_futures]

    coupling_futures = [extract_couplings(log, inputs=[]) for log in log_files]
    coupling_data = [f.result() for f in coupling_futures]

    combine_dir = os.path.abspath(f"{base_dir}/combined")
    combined_inp = create_combined_input(
        coupling_data=coupling_data,
        calc_params=calc_params,
        output_dir=combine_dir,
        inputs=[],
    ).result()

    combined_log = run_molcas(
        combined_inp, workdir_base, molcas_nprocs, inputs=[]
    ).result().strip()

    return {
        'name': name,
        'job_number': job_num,
        'n_roots': n_roots,
        'combined_dir': combine_dir,
        'combined_log': combined_log,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Po2 EFFE MS-CASPT2 workflow — 6 spin/symmetry blocks, C2 symmetry'
    )
    parser.add_argument('config', help='Path to config_Po2.yml')
    parser.add_argument(
        '--full-rasscf', action='store_true',
        help='Run full RASSCF optimization per spin group before CASPT2 '
             '(orbital handoff: singlet → triplet → quintet). '
             'Default: use pre-converged autoCAS orbitals directly (CIONLY).'
    )
    args = parser.parse_args()

    print(f"Loading config: {args.config}")
    config = load_config(args.config)

    print("Setting up Parsl (SLURM provider)...")
    setup_parsl(config['cluster'])

    xyz_file = config['molecule']['xyz_file']
    print(f"Loading geometry: {xyz_file}")
    with open(xyz_file, 'r') as f:
        xyz_content = f.read()

    workdir_base = config['cluster']['workdir_base']
    molcas_nprocs = config['cluster'].get('molcas_nprocs', 1)
    basis = config['calculations'][0]['basis']

    n_blocks = len(config['calculations'])
    total_roots = sum(c['n_roots'] for c in config['calculations'])
    print(f"\n{'='*60}")
    print(f"Po2 EFFE workflow: {n_blocks} blocks, {total_roots} root jobs total")
    print(f"Mode: {'full RASSCF + orbital handoff' if args.full_rasscf else 'CIONLY (pre-converged orbs)'}")
    print(f"{'='*60}")

    block_results = []

    if args.full_rasscf:
        # Group blocks by spin multiplicity (preserving config order)
        from itertools import groupby
        spin_groups = []
        for spin_val, group in groupby(config['calculations'], key=lambda b: b['spin']):
            spin_groups.append(list(group))

        # Use the first block's rasorb_path as the global starting orbital
        start_rasorb = os.path.abspath(config['calculations'][0]['rasorb_path'])
        print(f"Starting orbital: {start_rasorb}")

        for spin_group in spin_groups:
            results, start_rasorb = process_spin_group_with_rasscf(
                blocks=spin_group,
                xyz_content=xyz_content,
                workdir_base=workdir_base,
                molcas_nprocs=molcas_nprocs,
                start_rasorb=start_rasorb,
            )
            block_results.extend(results)
    else:
        # Original CIONLY mode: process each block sequentially
        for spin_config in config['calculations']:
            result = process_block(
                spin_config, xyz_content, workdir_base, molcas_nprocs
            )
            block_results.append(result)

    print(f"\n{'='*60}")
    print(f"All {n_blocks} blocks complete. Building final SO-RASSI...")
    print(f"{'='*60}")

    rassi_dir = "./final_rassi"
    rassi_inp = create_final_rassi_input(
        block_results=block_results,
        output_dir=rassi_dir,
        xyz_content=xyz_content,
        basis=basis,
        inputs=[],
    ).result()
    print(f"Final RASSI input: {rassi_inp}")

    print("Running final SO-RASSI (Spin Orbit, no Ejob)...")
    rassi_log = run_molcas(rassi_inp, workdir_base, molcas_nprocs, inputs=[]).result().strip()
    print(f"Final SO-RASSI done: {rassi_log}")

    print(f"\n{'='*60}")
    print("PO2 EFFE WORKFLOW COMPLETED SUCCESSFULLY")
    print(f"{'='*60}")
    for br in sorted(block_results, key=lambda x: x['job_number']):
        print(f"  JOB{br['job_number']:03d}: {br['name']} ({br['n_roots']} roots)")
    print(f"  Final RASSI log: {rassi_log}")
    print(f"{'='*60}")

    import subprocess
    result = subprocess.run(
        ['grep', '-m', '10', 'SO-RASSI State', rassi_log],
        capture_output=True, text=True
    )
    if result.stdout:
        print("\nSO-RASSI energies (first 10 states):")
        print(result.stdout)

    parsl.clear()


if __name__ == "__main__":
    main()
