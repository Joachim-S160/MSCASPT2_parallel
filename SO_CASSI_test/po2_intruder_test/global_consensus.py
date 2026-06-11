#!/usr/bin/env python3
"""
Po2 intruder-state test — global KILL_INTRUDERS consensus across all geometries.

After all per-geometry CASPT2 workers finish, this script:
  1. Reads every worker_{spin}_k*_b*.log from run_00/ … run_NN/ for each spin block.
  2. For each root index, finds the MINIMUM reference weight across all geometries.
  3. Finds N_surviving_global = the largest root index where ALL geometries pass the threshold.
     (First failure terminates — all higher roots are also dropped.)
  4. Writes global_consensus.json  — consumed by run_assemble.pbs and run_rassi.pbs.
  5. Writes refwt_all.dat           — human-readable long-format table for QA.

The globally consistent N_surviving is critical for dissociation curves: all RASSI jobs
must use the same state count so the SO energies are comparable across geometries.

Usage:
  cd po2_intruder_test/
  python3 global_consensus.py --n-geom 8 --kill-intruders 0.80
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Reuse extraction helpers from collect_heff.py (symlinked in same dir)
_SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(_SCRIPT_DIR))
from collect_heff import (
    extract_ref_weights,
    find_truncation_point,
    load_spin_logs,
)

# ---------------------------------------------------------------------------
# Spin configuration
# ---------------------------------------------------------------------------

_SPIN_LABEL = {"quintet": "Q", "triplet": "T", "singlet": "S"}


def _flag(w: float) -> str:
    if w >= 0.90:
        return "good"
    if w >= 0.80:
        return "okayish"
    if w >= 0.70:
        return "warning"
    return "severe"


# ---------------------------------------------------------------------------
# Core consensus computation
# ---------------------------------------------------------------------------

def compute_global_consensus(
    workdir: Path,
    n_geom: int,
    spin: str,
    n_roots: int,
    threshold: float,
) -> tuple[int, list[float], dict[str, list[float]]]:
    """
    Return (n_surviving_global, min_weights_per_root, per_geom_ref_weights).

    n_surviving_global = largest root index (1-based) where every geometry passes
    the threshold. First root below threshold terminates the count.
    """
    min_weights = [1.0] * n_roots
    per_geom: dict[str, list[float]] = {}

    for geom_idx in range(n_geom):
        geom_dir = workdir / f"run_{geom_idx:02d}"
        try:
            combined_log, _ = load_spin_logs(spin, geom_dir)
        except SystemExit:
            print(f"  ERROR: missing worker logs for {geom_dir.name} / {spin}",
                  file=sys.stderr)
            sys.exit(1)

        ref_weights = extract_ref_weights(combined_log, n_roots)
        per_geom[f"run_{geom_idx:02d}"] = ref_weights
        for i, w in enumerate(ref_weights):
            if w < min_weights[i]:
                min_weights[i] = w

    trunc_at = find_truncation_point(min_weights, threshold)
    if trunc_at == -1:
        n_surviving = n_roots
    else:
        n_surviving = trunc_at - 1

    return n_surviving, min_weights, per_geom


# ---------------------------------------------------------------------------
# refwt_all.dat writer
# ---------------------------------------------------------------------------

def write_refwt_all(
    output_path: Path,
    workdir: Path,
    n_geom: int,
    geom_r: list[float],
    spins: list[str],
    spin_roots: dict[str, int],
    per_geom_all: dict[str, dict[str, list[float]]],
) -> None:
    """
    Write refwt_all.dat in the same format as production_benchmark_revised
    07_extract_energies.sh:

      # Distance(A)  Spin  Root  RefWeight  Flag
      1.950000        Q     1    0.93971    good
    """
    lines = [
        "# CASPT2 reference weights for ALL roots of ALL spin manifolds (long format)",
        "# good>=0.90  okayish=0.80-0.90  warning=0.70-0.80  severe<0.70",
        "# Distance(A)  Spin  Root  RefWeight  Flag",
    ]
    for geom_idx in range(n_geom):
        key = f"run_{geom_idx:02d}"
        r = geom_r[geom_idx]
        for spin in spins:
            label = _SPIN_LABEL[spin]
            weights = per_geom_all[spin][key]
            for root_idx, w in enumerate(weights):
                root = root_idx + 1
                flag = _flag(w)
                lines.append(f"{r:<12.6f}  {label}  {root:>3d}  {w:.5f}  {flag}")

    output_path.write_text("\n".join(lines) + "\n")
    print(f"  Written: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute global KILL_INTRUDERS consensus across all geometries. "
                    "Writes global_consensus.json and refwt_all.dat."
    )
    parser.add_argument("--n-geom",          type=int, required=True,
                        help="Number of geometries (= number of run_NN/ directories)")
    parser.add_argument("--kill-intruders",  type=float, default=0.80,
                        help="Reference weight threshold (default: 0.80)")
    parser.add_argument("--ciroot-q",        type=int, required=True)
    parser.add_argument("--ciroot-t",        type=int, required=True)
    parser.add_argument("--ciroot-s",        type=int, required=True)
    parser.add_argument("--geom-r",          type=float, nargs="+", required=True,
                        help="R(Po-Po) in Å for each geometry (space-separated, same order "
                             "as run_00 … run_NN)")
    parser.add_argument("--workdir",         default=".",
                        help="Directory containing run_NN/ subdirs (default: .)")
    args = parser.parse_args()

    workdir   = Path(args.workdir).resolve()
    n_geom    = args.n_geom
    threshold = args.kill_intruders
    geom_r    = args.geom_r

    if len(geom_r) != n_geom:
        print(f"ERROR: --geom-r has {len(geom_r)} values but --n-geom={n_geom}",
              file=sys.stderr)
        sys.exit(1)

    spin_roots = {
        "quintet": args.ciroot_q,
        "triplet": args.ciroot_t,
        "singlet": args.ciroot_s,
    }
    spins = ["quintet", "triplet", "singlet"]

    print(f"\n{'='*60}")
    print(f"  Global KILL_INTRUDERS Consensus")
    print(f"  threshold={threshold}  n_geom={n_geom}")
    print(f"  CIROOT: Q={args.ciroot_q}  T={args.ciroot_t}  S={args.ciroot_s}")
    print(f"{'='*60}\n")

    consensus_result: dict = {"threshold": threshold}
    per_geom_all: dict[str, dict[str, list[float]]] = {}

    for spin in spins:
        n_roots = spin_roots[spin]
        label = _SPIN_LABEL[spin]
        print(f"--- {spin} ({label}, {n_roots} roots) ---")

        n_surviving, min_weights, per_geom = compute_global_consensus(
            workdir, n_geom, spin, n_roots, threshold
        )
        per_geom_all[spin] = per_geom

        # Per-geometry n_surviving for diagnostic
        per_geom_n_surv: dict[str, int] = {}
        for key, weights in per_geom.items():
            t = find_truncation_point(weights, threshold)
            per_geom_n_surv[key] = (t - 1) if t != -1 else n_roots

        first_fail = (n_surviving + 1) if n_surviving < n_roots else -1
        if first_fail == -1:
            print(f"  All {n_roots} roots pass at every geometry.")
        else:
            # Find which geometry first caused the failure
            worst_geom = min(
                per_geom.keys(),
                key=lambda k: per_geom[k][first_fail - 1]
            )
            worst_w = per_geom[worst_geom][first_fail - 1]
            print(f"  First failing root: {first_fail}  "
                  f"(worst geom: {worst_geom}, ref_weight={worst_w:.5f})")
        print(f"  N_surviving_global: {n_surviving} / {n_roots}")
        print(f"  Per-geom n_surviving: "
              + "  ".join(f"{k}={v}" for k, v in per_geom_n_surv.items()))
        print()

        consensus_result[f"n_surviving_{spin}"] = n_surviving
        consensus_result[f"per_spin_{spin}"] = {
            "n_roots": n_roots,
            "n_surviving": n_surviving,
            "first_failing_root": first_fail,
            "min_ref_weights": [round(w, 6) for w in min_weights],
            "per_geom_n_surviving": per_geom_n_surv,
        }

    # Write global_consensus.json
    json_path = workdir / "global_consensus.json"
    json_path.write_text(json.dumps(consensus_result, indent=2))
    print(f"Written: {json_path}")

    # Write refwt_all.dat
    refwt_path = workdir / "refwt_all.dat"
    write_refwt_all(refwt_path, workdir, n_geom, geom_r, spins, spin_roots, per_geom_all)

    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"  N_surviving: Q={consensus_result['n_surviving_quintet']}  "
          f"T={consensus_result['n_surviving_triplet']}  "
          f"S={consensus_result['n_surviving_singlet']}")
    print(f"  (RASSI will use: Nr of JobIphs = 3 "
          f"{consensus_result['n_surviving_quintet']} "
          f"{consensus_result['n_surviving_triplet']} "
          f"{consensus_result['n_surviving_singlet']})")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
