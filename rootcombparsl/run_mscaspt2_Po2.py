#!/usr/bin/env python3
"""
Po2 MS-CASPT2 EFFE workflow — NoSym SO-CASSI.

Pipelined --full-rasscf mode (default for production):
  Singlet RASSCF → submit singlet CASPT2 jobs + Triplet RASSCF
                 → submit triplet CASPT2 jobs + Quintet RASSCF
                 → submit quintet CASPT2 jobs
                 → collect all + EFFE combine + SO-RASSI

CASPT2 root jobs for each spin block are submitted to SLURM immediately
after that block's RASSCF completes, and run in parallel while the next
spin's RASSCF runs. Input files are created locally (fast) — not in SLURM.

CIONLY mode (--no-rasscf): uses pre-converged autoCAS orbitals directly,
skips RASSCF re-optimisation. Useful for testing.

Config: config_Po2.yml
  3 blocks: quintet (70 roots, JOB001), triplet (90 roots, JOB002),
            singlet (28 roots, JOB003). NoSym, CAS(12,8).
"""

import parsl
from parsl import bash_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider
from parsl.launchers import SimpleLauncher
import os
import re
import yaml
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import numpy as np
import argparse


SEWARD_FILES = ["RunFile", "OneInt", "SymInfo", "NqGrid",
                "ChVec1", "ChSel1", "ChRed", "ChDiag", "ChMap", "ChRst"]


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def load_config(config_file: str) -> Dict[str, Any]:
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)


def setup_parsl(cluster_config: Dict[str, Any]) -> None:
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
                    launcher=SimpleLauncher(),
                ),
            )
        ],
        strategy='htex_auto_scale',
        retries=1,
    )
    parsl.load(config)


# ---------------------------------------------------------------------------
# Input file creation — runs locally in master process (no Parsl overhead)
# ---------------------------------------------------------------------------

def create_root_input(root_idx: int, calc_params: Dict[str, Any], output_dir: str,
                      seward_dir: str = "") -> str:
    """Write RASSCF(CIONLY)+CASPT2(only=N) input for one root. Returns input file path."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    rasorb_path = os.path.abspath(calc_params['rasorb_path'])

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    if seward_dir:
        seward_section = "\n".join(
            f">>COPY {seward_dir}/seward.{f} $WorkDir/$Project.{f}"
            for f in SEWARD_FILES
        ) + "\n"
    else:
        seward_section = "&SEWARD\nCholesky\n"

    content = f"""&GATEWAY
Title
Po2 Root {root_idx}/{calc_params['n_roots']} spin={calc_params['spin']} sym={calc_params['symmetry']}
COORD = {xyz_file}
GROUP = NoSym
BASIS = {calc_params['basis']}
{seward_section}&RASSCF
File = {rasorb_path}
CIONLY
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {calc_params['n_roots']} {calc_params['n_roots']} 1
nActEl = {calc_params['nactel']} 0 0
Inactive = {calc_params['inactive']}
Ras2 = {calc_params['ras2']}
THRS = 1.0e-08 1.0e-04 1.0e-04
Levshft = 0.1
Iteration = 200 50
CIMX = 200
SDAV = 500
ORBAppear = COMPACT
&CASPT2
MAXITER = 300
Frozen = {calc_params['inactive']}
Multistate = all
Imaginary Shift = {calc_params['imaginary']}
only = {root_idx}
"""
    input_file = f"{output_dir}/root{root_idx}.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    return input_file


def create_caspt2_only_input(root_idx: int, calc_params: Dict[str, Any], output_dir: str,
                             seward_dir: str, jobmix_path: str) -> str:
    """Write CASPT2(only=N) input using injected JobIph — no RASSCF. Returns input file path."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    jobmix_path = os.path.abspath(jobmix_path)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    if seward_dir:
        seward_section = "\n".join(
            f">>COPY {seward_dir}/seward.{f} $WorkDir/$Project.{f}"
            for f in SEWARD_FILES
        ) + "\n"
    else:
        seward_section = "&SEWARD\nCholesky\n"

    content = f"""&GATEWAY
Title
Po2 Root {root_idx}/{calc_params['n_roots']} spin={calc_params['spin']} sym={calc_params['symmetry']}
COORD = {xyz_file}
GROUP = NoSym
BASIS = {calc_params['basis']}
{seward_section}>>COPY {jobmix_path} $WorkDir/$Project.JobIph
&CASPT2
MAXITER = 300
Frozen = {calc_params['inactive']}
Multistate = all
Imaginary Shift = {calc_params['imaginary']}
only = {root_idx}
"""
    input_file = f"{output_dir}/root{root_idx}.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    return input_file


def create_rasscf_input(calc_params: Dict[str, Any], start_rasorb: str, output_dir: str,
                        seward_dir: str = "") -> str:
    """Write full RASSCF input (no CIONLY) for orbital optimisation. Returns input file path."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)
    rasorb_path = os.path.abspath(start_rasorb)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    name = calc_params['name']
    if seward_dir:
        seward_section = "\n".join(
            f">>COPY {seward_dir}/seward.{f} $WorkDir/$Project.{f}"
            for f in SEWARD_FILES
        ) + "\n"
        save_section = f">>COPY $Project.JobIph $CurrDir/{name}.JobIph\n"
    else:
        seward_section = "&SEWARD\nCholesky\n"
        save_section = (
            "\n".join(f">>COPY $Project.{f} $CurrDir/seward.{f}" for f in SEWARD_FILES)
            + f"\n>>COPY $Project.JobIph $CurrDir/{name}.JobIph\n"
        )

    content = f"""&GATEWAY
Title
Po2 RASSCF spin={calc_params['spin']} sym={calc_params['symmetry']} (full opt)
COORD = {xyz_file}
GROUP = NoSym
BASIS = {calc_params['basis']}
{seward_section}&RASSCF
File = {rasorb_path}
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {calc_params['n_roots']} {calc_params['n_roots']} 1
nActEl = {calc_params['nactel']} 0 0
Inactive = {calc_params['inactive']}
Ras2 = {calc_params['ras2']}
THRS = 1.0e-08 1.0e-04 1.0e-04
Levshft = 0.1
Iteration = 200 50
CIMX = 200
SDAV = 500
>>COPY $Project.RasOrb $CurrDir/$Project.RasOrb
{save_section}"""
    input_file = f"{output_dir}/rasscf_full.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    return input_file


def create_combined_input(
    coupling_data: List[Tuple[int, List[float]]],
    calc_params: Dict[str, Any],
    output_dir: str,
    seward_dir: str = "",
    jobmix_path: str = "",
) -> str:
    """Write CASPT2(EFFE) combined input. Returns input file path.

    When jobmix_path is given, injects JobIph and skips RASSCF CIONLY entirely.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)

    coupling_data = sorted(coupling_data, key=lambda x: x[0])
    n_roots = calc_params['n_roots']

    H_eff = np.zeros((n_roots, n_roots))
    for root_idx, couplings in coupling_data:
        if len(couplings) != n_roots:
            raise ValueError(f"Root {root_idx}: got {len(couplings)} couplings, expected {n_roots}")
        H_eff[root_idx - 1, :] = couplings

    effe_block = "\n".join(
        [str(n_roots)] +
        [" ".join(f"{H_eff[i, j]:.14E}" for j in range(n_roots)) for i in range(n_roots)]
    )

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(calc_params['xyz_content'])

    if seward_dir:
        seward_section = "\n".join(
            f">>COPY {seward_dir}/seward.{f} $WorkDir/$Project.{f}"
            for f in SEWARD_FILES
        ) + "\n"
    else:
        seward_section = "&SEWARD\nCholesky\n"

    job_label = f"JOB{calc_params['job_number']:03d}"

    if jobmix_path:
        jobmix_path = os.path.abspath(jobmix_path)
        rasscf_block = f">>COPY {jobmix_path} $WorkDir/$Project.JobIph\n"
    else:
        rasorb_path = os.path.abspath(calc_params['rasorb_path'])
        rasscf_block = f"""&RASSCF
File = {rasorb_path}
CIONLY
SPIN = {calc_params['spin']}
Symmetry = {calc_params['symmetry']}
CIROOT = {n_roots} {n_roots} 1
nActEl = {calc_params['nactel']} 0 0
Inactive = {calc_params['inactive']}
Ras2 = {calc_params['ras2']}
THRS = 1.0e-08 1.0e-04 1.0e-04
Levshft = 0.1
Iteration = 200 50
CIMX = 200
SDAV = 500
ORBAppear = COMPACT
"""

    content = f"""&GATEWAY
Title
Po2 Combined EFFE spin={calc_params['spin']} sym={calc_params['symmetry']} -> {job_label}
COORD = {xyz_file}
GROUP = NoSym
BASIS = {calc_params['basis']}
{seward_section}{rasscf_block}&CASPT2
MAXITER = 300
Frozen = {calc_params['inactive']}
Multistate = all
Imaginary Shift = {calc_params['imaginary']}
EFFE
{effe_block}
>>COPY $Project.JobMix $CurrDir/{job_label}
"""
    input_file = f"{output_dir}/combined_effe.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    return input_file


def create_final_rassi_input(
    block_results: List[Dict[str, Any]],
    output_dir: str,
    xyz_content: str,
    basis: str,
    seward_dir: str = "",
) -> str:
    """Write final SO-RASSI input combining all JOBxxx files. Returns input file path."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_dir = os.path.abspath(output_dir)

    block_results = sorted(block_results, key=lambda x: x['job_number'])
    n_blocks = len(block_results)

    xyz_file = f"{output_dir}/molecule.xyz"
    with open(xyz_file, 'w') as f:
        f.write(xyz_content)

    copy_block = "\n".join(
        f">>COPY {blk['combined_dir']}/JOB{blk['job_number']:03d} JOB{blk['job_number']:03d}"
        for blk in block_results
    )

    if seward_dir:
        seward_section = "\n".join(
            f">>COPY {seward_dir}/seward.{f} $WorkDir/$Project.{f}"
            for f in SEWARD_FILES
        ) + "\n"
    else:
        seward_section = "&SEWARD\nCholesky\n"

    content = f"""&GATEWAY
Title
Po2 Final SO-RASSI {n_blocks} blocks
COORD = {xyz_file}
GROUP = NoSym
BASIS = {basis}
{seward_section}{copy_block}
&RASSI
Nr of JobIphs = {n_blocks} all
Spin Orbit
Ejob
Omega
End of Input
"""
    input_file = f"{output_dir}/final_rassi.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    return input_file


# ---------------------------------------------------------------------------
# SLURM jobs — only actual MOLCAS runs go through Parsl
# ---------------------------------------------------------------------------

@bash_app
def run_molcas(input_file: str, workdir_base: str, nprocs: int, inputs=[]) -> str:
    """Submit one MOLCAS job to SLURM. Returns exit code (int) via .result()."""
    import os
    import time
    from pathlib import Path

    input_path = Path(input_file)
    output_dir = str(input_path.parent)
    job_name = input_path.stem
    workdir = f"{workdir_base}/parsl_{os.getpid()}_{int(time.time())}_{job_name}"

    return f"""
set -e
export MOLCAS_WORKDIR={workdir}
mkdir -p $MOLCAS_WORKDIR
cd {output_dir}
pymolcas -np {nprocs} {input_file} -f
"""


@bash_app
def run_rasscf_and_copy_orb(
    input_file: str, workdir_base: str, nprocs: int, out_rasorb: str, inputs=[]
) -> str:
    """Run full RASSCF and copy output .RasOrb to a stable path. Returns that path via stdout."""
    import os
    import time
    from pathlib import Path

    input_path = Path(input_file)
    output_dir = str(input_path.parent)
    job_name = input_path.stem
    workdir = f"{workdir_base}/parsl_{os.getpid()}_{int(time.time())}_{job_name}"

    return f"""
set -e
export MOLCAS_WORKDIR={workdir}
mkdir -p $MOLCAS_WORKDIR
cd {output_dir}
pymolcas -np {nprocs} {input_file} -f
if [ ! -f {output_dir}/{job_name}.RasOrb ]; then
    echo "ERROR: {output_dir}/{job_name}.RasOrb not found after RASSCF"
    ls -la {output_dir}/
    exit 1
fi
cp {output_dir}/{job_name}.RasOrb {out_rasorb}
echo {out_rasorb}
"""


# ---------------------------------------------------------------------------
# Local helpers (run in master process)
# ---------------------------------------------------------------------------

def extract_couplings(log_file: str) -> Tuple[int, List[float]]:
    """Parse 'Hamiltonian Effective Couplings' from a CASPT2(only=N) log."""
    root_match = re.search(r'root(\d+)', Path(log_file).stem)
    if not root_match:
        raise ValueError(f"Cannot extract root index from: {Path(log_file).stem}")
    root_idx = int(root_match.group(1))

    with open(log_file, 'r') as f:
        content = f.read()

    match = re.search(
        r'Hamiltonian Effective Couplings.*?\n((?:.*?<.*?\n)*)',
        content, re.DOTALL | re.IGNORECASE
    )
    if not match:
        raise ValueError(f"Could not find 'Hamiltonian Effective Couplings' in {log_file}")

    couplings = [float(line.split()[-1]) for line in match.group(1).split('\n') if '<' in line]
    return (root_idx, couplings)


def _log_is_complete(log_path: str) -> bool:
    """Return True if log exists and contains OpenMolcas completion marker."""
    if not os.path.isfile(log_path):
        return False
    try:
        with open(log_path, 'r') as f:
            return 'Happy landing' in f.read()
    except OSError:
        return False


def _make_calc_params(spin_config: Dict[str, Any], xyz_content: str,
                      rasorb_path: str) -> Dict[str, Any]:
    return {
        'name':       spin_config['name'],
        'rasorb_path': rasorb_path,
        'xyz_content': xyz_content,
        'n_roots':    spin_config['n_roots'],
        'spin':       spin_config['spin'],
        'symmetry':   spin_config['symmetry'],
        'inactive':   spin_config['inactive'],
        'ras2':       spin_config['ras2'],
        'nactel':     spin_config['nactel'],
        'basis':      spin_config.get('basis', 'ANO-RCC-VQZP'),
        'imaginary':  spin_config.get('imaginary', 0.25),
        'job_number': spin_config['job_number'],
    }


# ---------------------------------------------------------------------------
# Workflow building blocks
# ---------------------------------------------------------------------------

def _run_rasscf_for_block(
    spin_config: Dict[str, Any],
    start_rasorb: str,
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int,
    restart: bool = False,
    seward_dir: str = "",
) -> Tuple[str, str]:
    """Run full RASSCF for one spin block. Blocks until done.
    Returns (rasorb_path, jobiph_path)."""
    name = spin_config['name']
    rasscf_dir = os.path.abspath(f"{spin_config['output_dir']}/rasscf")
    out_rasorb = os.path.abspath(f"{rasscf_dir}/{name}.RasOrb")
    out_jobiph = os.path.abspath(f"{rasscf_dir}/{name}.JobIph")

    if restart and os.path.isfile(out_rasorb) and os.path.isfile(out_jobiph):
        print(f"  RASSCF [{name}] — skipping (RasOrb+JobIph exist)")
        return out_rasorb, out_jobiph

    calc_params = _make_calc_params(spin_config, xyz_content, start_rasorb)
    inp = create_rasscf_input(calc_params, start_rasorb, rasscf_dir, seward_dir=seward_dir)
    run_rasscf_and_copy_orb(
        inp, workdir_base, molcas_nprocs, out_rasorb, inputs=[]
    ).result()  # blocks; raises BashExitFailure on non-zero exit
    print(f"  RASSCF [{name}] complete -> {out_rasorb}")
    return out_rasorb, out_jobiph


def _launch_caspt2_jobs(
    spin_config: Dict[str, Any],
    rasorb_path: str,
    xyz_content: str,
    workdir_base: str,
    molcas_nprocs: int,
    restart: bool = False,
    seward_dir: str = "",
    jobmix_path: str = "",
) -> Tuple[Dict[str, Any], List, str]:
    """Create input files locally, then submit CASPT2 root jobs to SLURM (non-blocking).

    When jobmix_path is given, each root job injects the shared JobIph and runs only
    &CASPT2 (no RASSCF CIONLY). Otherwise falls back to RASSCF CIONLY + CASPT2.

    Returns (calc_params, mol_futures, log_files, combine_dir) for later collection.
    mol_futures[i] is None if restart=True and that root's log is already complete.
    """
    name = spin_config['name']
    n_roots = spin_config['n_roots']
    base_dir = spin_config['output_dir']
    calc_params = _make_calc_params(spin_config, xyz_content, rasorb_path)

    print(f"  Writing {n_roots} input files [{name}] locally...")
    if jobmix_path:
        input_files = [
            create_caspt2_only_input(r, calc_params, f"{base_dir}/root{r}", seward_dir, jobmix_path)
            for r in range(1, n_roots + 1)
        ]
    else:
        input_files = [
            create_root_input(r, calc_params, f"{base_dir}/root{r}", seward_dir=seward_dir)
            for r in range(1, n_roots + 1)
        ]

    log_files = [str(Path(inp).with_suffix('.log')) for inp in input_files]

    if restart:
        skipped = sum(1 for log in log_files if _log_is_complete(log))
        to_submit = n_roots - skipped
        print(f"  Restart: {skipped}/{n_roots} roots already complete [{name}], submitting {to_submit}")
    else:
        print(f"  Submitting {n_roots} CASPT2 jobs [{name}] to SLURM (running in background)...")

    mol_futures = [
        None if (restart and _log_is_complete(log_files[i]))
        else run_molcas(inp, workdir_base, molcas_nprocs, inputs=[])
        for i, inp in enumerate(input_files)
    ]

    return calc_params, mol_futures, log_files, os.path.abspath(f"{base_dir}/combined")


def _collect_and_effe(
    calc_params: Dict[str, Any],
    mol_futures: List,
    log_files: List[str],
    combine_dir: str,
    workdir_base: str,
    molcas_nprocs: int,
    restart: bool = False,
    seward_dir: str = "",
    jobmix_path: str = "",
) -> Dict[str, Any]:
    """Wait for CASPT2 root jobs, extract couplings, run EFFE combine."""
    name = calc_params['name']
    n_roots = calc_params['n_roots']
    job_num = calc_params['job_number']

    print(f"\n--- Waiting for {n_roots} CASPT2 jobs [{name}] ---")
    for f in mol_futures:
        if f is not None:
            f.result()  # blocks; raises BashExitFailure on non-zero exit
    print(f"  All jobs done [{name}]. Extracting couplings...")

    coupling_data = [extract_couplings(log) for log in log_files]
    print(f"  Couplings extracted [{name}]. Running EFFE combine...")

    job_file = os.path.join(combine_dir, f"JOB{job_num:03d}")
    combined_log = str(Path(combine_dir) / 'combined_effe.log')
    if restart and os.path.isfile(job_file):
        print(f"  EFFE [{name}] — skipping (JOB{job_num:03d} exists)")
        return {'name': name, 'job_number': job_num, 'n_roots': n_roots,
                'combined_dir': combine_dir, 'combined_log': combined_log}

    combined_inp = create_combined_input(coupling_data, calc_params, combine_dir,
                                         seward_dir=seward_dir, jobmix_path=jobmix_path)
    run_molcas(combined_inp, workdir_base, molcas_nprocs, inputs=[]).result()
    print(f"  EFFE done [{name}]: {combined_log}")

    return {'name': name, 'job_number': job_num, 'n_roots': n_roots,
            'combined_dir': combine_dir, 'combined_log': combined_log}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Po2 EFFE MS-CASPT2 workflow — NoSym, 3 spin blocks')
    parser.add_argument('config', help='Path to config_Po2.yml')
    parser.add_argument(
        '--no-rasscf', action='store_true',
        help='CIONLY mode: skip RASSCF re-optimisation, use autoCAS orbitals directly. '
             'Default: full RASSCF per spin block with pipelined CASPT2.'
    )
    parser.add_argument(
        '--restart', action='store_true',
        help='Resume a partially completed run: skip steps whose outputs already exist '
             '(RasOrb files, root CASPT2 logs with "Happy landing", JOB/SO-RASSI files).'
    )
    args = parser.parse_args()
    full_rasscf = not args.no_rasscf
    restart = args.restart

    config = load_config(args.config)
    setup_parsl(config['cluster'])

    with open(config['molecule']['xyz_file'], 'r') as f:
        xyz_content = f.read()

    workdir_base = config['cluster']['workdir_base']
    molcas_nprocs = config['cluster'].get('molcas_nprocs', 1)
    basis = config['calculations'][0]['basis']
    n_blocks = len(config['calculations'])
    total_roots = sum(c['n_roots'] for c in config['calculations'])

    print(f"\n{'='*60}")
    print(f"Po2 EFFE: {n_blocks} blocks, {total_roots} root jobs")
    print(f"Mode: {'full RASSCF + pipelined CASPT2' if full_rasscf else 'CIONLY'}"
          + (" [RESTART]" if restart else ""))
    print(f"{'='*60}\n")

    block_results = []
    seward_dir = ""  # set in full_rasscf branch; empty → &SEWARD in CIONLY mode

    if full_rasscf:
        # Process spin blocks in ascending order: singlet(1)→triplet(3)→quintet(5)
        # so RASSCF handoff uses the most closed-shell reference first.
        spin_order = sorted(config['calculations'], key=lambda b: b['spin'])
        current_rasorb = os.path.abspath(config['molecule']['rasorb_path'])
        print(f"Starting orbital: {current_rasorb}")

        # Singlet RASSCF saves all 10 Seward files to its rasscf dir; all subsequent
        # jobs inject them instead of rerunning &SEWARD (195× → 1× Seward total).
        singlet_config = min(spin_order, key=lambda b: b['spin'])
        seward_dir = os.path.abspath(f"{singlet_config['output_dir']}/rasscf")

        pending = []        # (calc_params, mol_futures, log_files, combine_dir)
        jobiph_paths = []   # per-block jobiph, aligned with pending

        for i, spin_config in enumerate(spin_order):
            name = spin_config['name']
            print(f"\n=== RASSCF [{name}] spin={spin_config['spin']}, {spin_config['n_roots']} roots ===")

            # For singlet (i=0): run &SEWARD to create files; for triplet/quintet: inject.
            rasscf_seward = (
                seward_dir
                if (i > 0 and os.path.isfile(os.path.join(seward_dir, 'seward.ChVec1')))
                else ""
            )

            # Run RASSCF — blocks. Previously submitted CASPT2 jobs run in parallel on SLURM.
            current_rasorb, jobiph = _run_rasscf_for_block(
                spin_config, current_rasorb, xyz_content, workdir_base, molcas_nprocs,
                restart=restart, seward_dir=rasscf_seward,
            )

            # Submit CASPT2 root jobs immediately — non-blocking.
            # They run on SLURM while the next spin's RASSCF executes.
            # By the time SLURM executes these, singlet RASSCF has already created seward+JobIph.
            pending.append(
                _launch_caspt2_jobs(
                    spin_config, current_rasorb, xyz_content, workdir_base, molcas_nprocs,
                    restart=restart, seward_dir=seward_dir, jobmix_path=jobiph,
                )
            )
            jobiph_paths.append(jobiph)

        print(f"\n{'='*60}")
        print("All RASSCF complete. Collecting CASPT2 + assembling EFFE...")
        print(f"{'='*60}")

        for (calc_params, mol_futures, log_files, combine_dir), jobiph in zip(pending, jobiph_paths):
            block_results.append(
                _collect_and_effe(calc_params, mol_futures, log_files, combine_dir, workdir_base, molcas_nprocs,
                                  restart=restart, seward_dir=seward_dir, jobmix_path=jobiph)
            )

    else:
        # CIONLY: process each block sequentially with pre-converged autoCAS orbitals
        rasorb = os.path.abspath(config['molecule']['rasorb_path'])
        for spin_config in config['calculations']:
            name = spin_config['name']
            n_roots = spin_config['n_roots']
            job_num = spin_config['job_number']
            calc_params = _make_calc_params(spin_config, xyz_content, rasorb)
            base_dir = spin_config['output_dir']

            print(f"\n--- Block [{name}] spin={spin_config['spin']}, {n_roots} roots ---")
            input_files = [
                create_root_input(r, calc_params, f"{base_dir}/root{r}")
                for r in range(1, n_roots + 1)
            ]
            log_files = [str(Path(inp).with_suffix('.log')) for inp in input_files]
            if restart:
                skipped = sum(1 for log in log_files if _log_is_complete(log))
                to_submit = n_roots - skipped
                print(f"  Restart: {skipped}/{n_roots} roots already complete [{name}], submitting {to_submit}")
            mol_futures = [
                None if (restart and _log_is_complete(log_files[i]))
                else run_molcas(inp, workdir_base, molcas_nprocs, inputs=[])
                for i, inp in enumerate(input_files)
            ]
            combine_dir = os.path.abspath(f"{base_dir}/combined")
            block_results.append(
                _collect_and_effe(calc_params, mol_futures, log_files, combine_dir, workdir_base, molcas_nprocs,
                                  restart=restart)
            )

    # Final SO-RASSI
    print(f"\n{'='*60}")
    print("Building final SO-RASSI...")
    rassi_dir = "./final_rassi"
    rassi_seward = seward_dir if (full_rasscf and os.path.isfile(os.path.join(seward_dir, 'seward.ChVec1'))) else ""
    rassi_inp = create_final_rassi_input(block_results, rassi_dir, xyz_content, basis, seward_dir=rassi_seward)
    rassi_log = str(Path(rassi_inp).with_suffix('.log'))
    if restart and _log_is_complete(rassi_log):
        print(f"Final SO-RASSI — skipping (log exists and complete)")
    else:
        run_molcas(rassi_inp, workdir_base, molcas_nprocs, inputs=[]).result()
        print(f"Final SO-RASSI done: {rassi_log}")

    print(f"\n{'='*60}")
    print("PO2 EFFE WORKFLOW COMPLETED")
    for br in sorted(block_results, key=lambda x: x['job_number']):
        print(f"  JOB{br['job_number']:03d}: {br['name']} ({br['n_roots']} roots)")
    print(f"  SO-RASSI log: {rassi_log}")
    print(f"{'='*60}")

    import subprocess
    grep = subprocess.run(['grep', '-m', '10', 'SO-RASSI State', rassi_log],
                         capture_output=True, text=True)
    if grep.stdout:
        print("\nSO-RASSI energies (first 10 states):")
        print(grep.stdout)

    parsl.clear()


if __name__ == "__main__":
    main()
