#!/usr/bin/env python3
"""
N2 CAS(6,6)/cc-pVDZ CASPT2 k-scaling prototype.

Measures how wall time and peak RSS scale with k (roots per CPU) to find
the sweet spot for the full Po2 MS-CASPT2 parallel benchmark on HPC.

Workflow:
  Phase 1 (once): GATEWAY + SEWARD + RASSCF → shared/ (11 files)
  Phase 2 (per k): run ceil(8/k) batches, each batch = k sequential ONLY= calls
    - measure wall time and peak RSS per batch via /usr/bin/time -v
  Report: table of k → per-batch wall, per-batch RSS, projected parallel wall

Usage:
  cd SO_CASSI_test/n2_scaling_test/
  python3 run_scaling_test.py [--k-values 1 2 4 8] [--n-roots 8] [--skip-setup]
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
BENCH_DIR = Path(__file__).parent.resolve()
MOLCAS_BUILD = Path("/home/joaschee/OpenMolcas/build")
PYMOLCAS     = MOLCAS_BUILD / "pymolcas"
MOLCAS_ENV   = Path("/home/joaschee/molcas_env/bin/activate")

SEWARD_FILES = ["RunFile", "OneInt", "SymInfo", "NqGrid",
                "ChVec1", "ChSel1", "ChRed", "ChDiag", "ChMap", "ChRst"]

# ---------------------------------------------------------------------------
# Input generation
# ---------------------------------------------------------------------------
SEWARD_INJECT = "\n".join(
    f">>COPY $CurrDir/shared/seward.{f} $WorkDir/$Project.{f}"
    for f in SEWARD_FILES
)

def make_batch_input(roots: list[int], frozen: int = 2,
                     shift: float = 0.25) -> str:
    """Return an OpenMolcas input string for a batch of ONLY= roots."""
    caspt2_blocks = "\n\n".join(f"""\
* --- Root {r} ---
&CASPT2
  MAXITER=  300
  Frozen=   {frozen}
  Multistate= all
  Imaginary Shift= {shift}
  Only= {r}
End of Input""" for r in roots)

    return f"""\
{SEWARD_INJECT}
>>COPY $CurrDir/shared/rasscf.JobIph $WorkDir/$Project.JobIph

{caspt2_blocks}
"""

# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def run_with_time(input_path: Path, workdir: Path, log_path: Path,
                  time_path: Path) -> tuple[float, float]:
    """
    Run pymolcas with /usr/bin/time -v.
    Returns (wall_s, peak_rss_mb).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["MOLCAS"] = str(MOLCAS_BUILD)
    env["MOLCAS_WORKDIR"] = str(workdir)
    env["MOLCAS_MEM"] = "2000"
    env["OMP_NUM_THREADS"] = "1"

    cmd = [
        "/usr/bin/time", "-v",
        str(PYMOLCAS), str(input_path),
    ]
    t0 = time.time()
    with log_path.open("w") as flog, time_path.open("w") as ftime:
        subprocess.run(
            cmd, stdout=flog, stderr=ftime,
            cwd=str(BENCH_DIR), env=env,
        )
    wall_s = time.time() - t0  # Python-level wall time (sub-second precision)

    # /usr/bin/time -v gives peak RSS — we only parse that from its output
    time_text = time_path.read_text(errors="replace")

    rss_mb = 0.0
    m = re.search(r"Maximum resident set size \(kbytes\):\s+(\d+)", time_text)
    if m:
        rss_mb = int(m.group(1)) / 1024
    else:
        # macOS /usr/bin/time uses bytes
        m = re.search(r"maximum resident set size\s+(\d+)", time_text, re.IGNORECASE)
        if m:
            rss_mb = int(m.group(1)) / (1024 * 1024)

    happy = "happy landing" in log_path.read_text(errors="replace").lower()
    if not happy:
        print(f"  WARNING: no 'Happy Landing' in {log_path.name}", file=sys.stderr)

    return wall_s, rss_mb


# ---------------------------------------------------------------------------
# Phase 1: setup run
# ---------------------------------------------------------------------------
def run_phase1(skip: bool) -> None:
    shared = BENCH_DIR / "shared"
    if skip and shared.exists() and (shared / "rasscf.JobIph").exists():
        print("Phase 1: skipping (shared/ already exists)")
        return

    print("Phase 1: GATEWAY + SEWARD + RASSCF ...")
    shared.mkdir(exist_ok=True)
    workdir = BENCH_DIR / "mw_setup"
    shutil.rmtree(workdir, ignore_errors=True)
    workdir.mkdir()

    wall, rss = run_with_time(
        input_path=BENCH_DIR / "phase1_setup.input",
        workdir=workdir,
        log_path=BENCH_DIR / "phase1.log",
        time_path=BENCH_DIR / "phase1_time.txt",
    )
    files = list(shared.iterdir())
    print(f"  Done: {wall:.1f} s  RSS={rss:.1f} MB  ({len(files)} files in shared/)")
    if not (shared / "rasscf.JobIph").exists():
        print("ERROR: rasscf.JobIph not found in shared/ — phase1 failed. Check phase1.log.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Phase 2: k-scaling
# ---------------------------------------------------------------------------
def run_k_scaling(k_values: list[int], n_roots: int) -> list[dict]:
    results = []

    for k in k_values:
        # Split n_roots into batches of size k
        batches = [
            list(range(i + 1, min(i + k, n_roots) + 1))
            for i in range(0, n_roots, k)
        ]
        n_batches = len(batches)
        print(f"\nk={k}: {n_batches} batch(es) × up to {k} root(s) each")

        batch_walls = []
        batch_rss   = []

        for b_idx, roots in enumerate(batches):
            input_text = make_batch_input(roots)
            input_path = BENCH_DIR / f"batch_k{k}_b{b_idx}.input"
            input_path.write_text(input_text)

            workdir  = BENCH_DIR / f"mw_k{k}_b{b_idx}"
            log_path = BENCH_DIR / f"log_k{k}_b{b_idx}.log"
            time_path = BENCH_DIR / f"time_k{k}_b{b_idx}.txt"

            shutil.rmtree(workdir, ignore_errors=True)
            workdir.mkdir()

            wall, rss = run_with_time(input_path, workdir, log_path, time_path)
            batch_walls.append(wall)
            batch_rss.append(rss)
            print(f"  batch {b_idx+1}/{n_batches} roots={roots}: "
                  f"wall={wall:.2f}s  RSS={rss:.1f}MB")

        # Projected parallel wall time = longest single batch (if all batches ran in parallel)
        parallel_wall = max(batch_walls)
        mean_rss      = sum(batch_rss) / len(batch_rss)
        max_rss       = max(batch_rss)

        results.append({
            "k": k,
            "n_batches": n_batches,
            "batch_walls_s": batch_walls,
            "batch_rss_mb": batch_rss,
            "parallel_wall_s": parallel_wall,    # bottleneck batch = total parallel time
            "mean_rss_mb": mean_rss,
            "max_rss_mb": max_rss,
        })

    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def print_report(results: list[dict], n_roots: int) -> None:
    print("\n" + "=" * 72)
    print(f"  N2 CAS(6,6)/cc-pVDZ — CASPT2 k-scaling (n_roots={n_roots})")
    print("=" * 72)
    print(f"  {'k':>4}  {'n_batches':>9}  {'wall/batch(s)':>13}  "
          f"{'proj_parallel(s)':>16}  {'mean_RSS(MB)':>12}  {'max_RSS(MB)':>11}")
    print("  " + "-" * 70)
    for r in results:
        print(f"  {r['k']:>4}  {r['n_batches']:>9}  "
              f"{r['parallel_wall_s']:>13.2f}  "
              f"{r['parallel_wall_s']:>16.2f}  "
              f"{r['mean_rss_mb']:>12.1f}  "
              f"{r['max_rss_mb']:>11.1f}")
    print("=" * 72)
    print()
    print("  Interpretation:")
    print("  - wall/batch = time for the bottleneck batch = total HPC wall if all batches run in parallel")
    print("  - RSS = peak memory per CPU (batch); total HPC memory = n_batches × max_RSS")
    print("  - Sweet spot: k where parallel wall is minimized and n_batches×RSS fits in node memory")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="N2 CASPT2 k-scaling prototype")
    parser.add_argument("--k-values", type=int, nargs="+", default=[1, 2, 4, 8],
                        help="Roots-per-batch values to test (default: 1 2 4 8)")
    parser.add_argument("--n-roots", type=int, default=8,
                        help="Total number of CASPT2 roots (default: 8)")
    parser.add_argument("--skip-setup", action="store_true",
                        help="Skip phase 1 if shared/ already exists")
    args = parser.parse_args()

    # Source molcas_env in subprocess env — already done at shell level if script
    # is invoked via: source /home/joaschee/molcas_env/bin/activate && python3 ...
    # but we set MOLCAS explicitly so pymolcas works regardless.
    os.environ["MOLCAS"] = str(MOLCAS_BUILD)

    run_phase1(skip=args.skip_setup)
    results = run_k_scaling(args.k_values, args.n_roots)
    print_report(results, args.n_roots)


if __name__ == "__main__":
    main()
