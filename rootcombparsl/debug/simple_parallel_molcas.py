#!/usr/bin/env python3
"""
Simple script to run 3 MOLCAS jobs in parallel using Parsl
This is a minimal test to verify the setup works
"""
import parsl
from parsl import bash_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider
import os
import re

# Setup Parsl configuration
config = Config(
    executors=[
        HighThroughputExecutor(
            label="molcas_executor",
            cores_per_worker=1,
            provider=SlurmProvider(
                partition="cpu_milan_rhel9",
                account="2025_026",
                nodes_per_block=1,
                init_blocks=0,
                max_blocks=3,  # Only 3 parallel jobs
                walltime="00:30:00",
                scheduler_options="#SBATCH --clusters=dodrio\n",
                worker_init="""
module purge
module load autoCAS/3.0.0-iomkl-2023a
export MOLCAS_PRINT=verbose
export MOLCAS_MEM=4000
export MOLCAS_NPROCS=1
export OMP_NUM_THREADS=1
ulimit -s unlimited

echo "=== Worker Environment ==="
echo "Hostname: $(hostname)"
echo "Date: $(date)"
which pymolcas
echo "=========================="
                """,
            ),
        )
    ],
    strategy='htex_auto_scale',
)

parsl.load(config)


@bash_app
def run_molcas_job(input_file: str, workdir_base: str) -> str:
    """Run a single MOLCAS job"""
    return f"""
    set -x
    
    echo "=== Starting MOLCAS job ==="
    echo "Input file: {input_file}"
    echo "Current directory: $(pwd)"
    echo "Hostname: $(hostname)"
    
    # Check if pymolcas is available
    which pymolcas || (echo "ERROR: pymolcas not found" && exit 1)
    
    # Get the directory and filename
    INPUT_DIR=$(dirname {input_file})
    INPUT_NAME=$(basename {input_file} .inp)
    
    # Create unique workdir
    WORKDIR="{workdir_base}/test_${{INPUT_NAME}}_$_$(date +%s)"
    export MOLCAS_WORKDIR=$WORKDIR
    mkdir -p $MOLCAS_WORKDIR
    
    echo "Working directory: $MOLCAS_WORKDIR"
    
    # Check if input file exists
    if [ ! -f "{input_file}" ]; then
        echo "ERROR: Input file {input_file} not found"
        exit 1
    fi
    
    # Change to input directory
    cd $INPUT_DIR
    echo "Changed to: $(pwd)"
    
    # Run pymolcas
    echo "Running pymolcas..."
    pymolcas -np 1 $(basename {input_file}) -f
    
    # Check for log file and "Happy landing"
    if [ ! -f "${{INPUT_NAME}}.log" ]; then
        echo "ERROR: Log file ${{INPUT_NAME}}.log not created"
        echo "Files in directory:"
        ls -la
        exit 1
    fi
    
    # Check for successful completion
    if grep -q "Happy landing" ${{INPUT_NAME}}.log; then
        echo "✓ MOLCAS completed successfully (Happy landing found)"
        # Output ONLY the log file path on the last line (this is what gets returned)
        echo "$(pwd)/${{INPUT_NAME}}.log"
        exit 0
    else
        echo "✗ MOLCAS may have failed (no Happy landing found)"
        echo "Last 30 lines of log:"
        tail -n 30 ${{INPUT_NAME}}.log
        exit 1
    fi
    """


def prepare_input_from_template(template_file: str, root_num: int, output_dir: str) -> str:
    """
    Create input file from template by replacing only=X with the root number
    
    Args:
        template_file: Path to template .inp file
        root_num: Root number (1, 2, 3, ...)
        output_dir: Directory to create the input file
    
    Returns:
        Absolute path to created input file
    """
    # Read template
    with open(template_file, 'r') as f:
        content = f.read()
    
    # Replace only=X with the correct root number
    # This regex finds "only=<number>" and replaces it
    content = re.sub(r'only\s*=\s*\d+', f'only={root_num}', content, flags=re.IGNORECASE)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Write input file
    input_file = f"{output_dir}/root{root_num}.inp"
    with open(input_file, 'w') as f:
        f.write(content)
    
    return os.path.abspath(input_file)


def main():
    print("=" * 60)
    print("Simple Parallel MOLCAS Test")
    print("=" * 60)
    
    # Configuration - UPDATE THESE
    template_file = "template_root.inp"  # Your template input file
    workdir_base = "/dodrio/scratch/projects/starting_2025_097/molcas_workdirs"
    n_roots = 3  # Number of roots to run
    
    # Check if template exists
    if not os.path.exists(template_file):
        print(f"ERROR: Template file not found: {template_file}")
        print("Please create your template_root.inp file first.")
        return
    
    print(f"Template file: {os.path.abspath(template_file)}")
    print(f"Workdir base: {workdir_base}")
    print(f"Number of roots: {n_roots}")
    print()
    
    # Create input files from template
    print(f"Creating {n_roots} input files from template...")
    input_files = []
    for i in range(1, n_roots + 1):
        output_dir = f"test_root{i}"
        inp = prepare_input_from_template(template_file, i, output_dir)
        input_files.append(inp)
        print(f"  Root {i}: {inp}")
    
    print()
    print(f"Submitting {n_roots} MOLCAS jobs in parallel...")
    print("(This will submit worker jobs to the cluster)")
    
    # Submit all jobs
    futures = []
    for inp in input_files:
        future = run_molcas_job(inp, workdir_base)
        futures.append(future)
    
    # Wait for results (just to ensure jobs complete, we don't need the return values)
    print()
    print("Waiting for jobs to complete...")
    print("(Check 'squeue -u $USER' to see worker jobs)")
    print()
    
    results = []
    for i, future in enumerate(futures, 1):
        try:
            print(f"Root {i}: Waiting...")
            # Just call result() to wait for completion, ignore the return value
            future.result()
            
            # Log file is always in a predictable location
            log_file = f"test_root{i}/root{i}.log"
            
            if not os.path.exists(log_file):
                raise Exception(f"Log file not found: {log_file}")
            
            # Extract couplings from log file
            with open(log_file, 'r') as f:
                log_content = f.read()
            
            # Find Hamiltonian Effective Couplings section
            match = re.search(r'Hamiltonian Effective Couplings.*?\n((?:.*?<.*?\n)*)', 
                            log_content, re.DOTALL | re.IGNORECASE)
            
            if not match:
                raise Exception("Could not find Effective Couplings in log file")
            
            # Extract coupling values
            coupling_lines = match.group(1)
            couplings = []
            for line in coupling_lines.split('\n'):
                if '<' in line:
                    # Extract the last column (coupling value)
                    value = float(line.split()[-1])
                    couplings.append(value)
            
            print(f"Root {i}: SUCCESS ✓")
            print(f"  Log: {os.path.abspath(log_file)}")
            print(f"  Couplings: {couplings}")
            results.append(("SUCCESS", log_file, couplings))
            
        except Exception as e:
            print(f"Root {i}: FAILED ✗")
            print(f"  Error: {e}")
            results.append(("FAILED", None, None))
    
    # Summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    successes = sum(1 for r in results if r[0] == "SUCCESS")
    print(f"Successful: {successes}/{n_roots}")
    print(f"Failed: {n_roots-successes}/{n_roots}")
    
    if successes == n_roots:
        print()
        print("✓ All jobs completed successfully!")
        print()
        print("Extracted Couplings:")
        print("-" * 60)
        for i, r in enumerate(results, 1):
            if r[0] == "SUCCESS":
                couplings = r[2]
                print(f"Root {i}: {couplings}")
        print("-" * 60)
    else:
        print()
        print("✗ Some jobs failed. Check the error messages above.")
    
    print("=" * 60)
    
    # Cleanup Parsl
    parsl.clear()


if __name__ == "__main__":
    main()
