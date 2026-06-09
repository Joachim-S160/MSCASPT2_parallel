#!/usr/bin/env python3
"""
Compare MS-CASPT2 state energies between Run A (Multistate=all) and
Run B (EFFE assembled).  Parses both logs and checks that all N state
energies match to within the specified threshold (default 1e-7 Ha).

Tries three energy patterns in order (first that yields N values wins):
  1. RASSI State   N  Total energy:   E      (non-SO RASSI output)
  2. SO-RASSI State N  Total energy:  E      (if SO-RASSI was used)
  3. MS-CASPT2 Root N  Total energy = E      (CASPT2-level eigenvalues)

Prints a per-root comparison table and a PASS / FAIL verdict.
"""

import argparse
import re
import sys


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_PATTERNS = [
    ("RASSI-state",   re.compile(r"(?<!SO-)RASSI State\s+(\d+)\s+Total energy:\s*([-\d.]+)", re.IGNORECASE)),
    ("SO-RASSI-state", re.compile(r"SO-RASSI State\s+(\d+)\s+Total energy:\s*([-\d.]+)", re.IGNORECASE)),
    ("MS-CASPT2-root", re.compile(r"MS-CASPT2 Root\s+(\d+)\s+Total energy\s*[=:]\s*([-\d.]+)", re.IGNORECASE)),
]


def parse_energies(log_text: str, n_roots: int) -> tuple[str, list[float]]:
    """Return (pattern_name, [energy_1, ..., energy_n]) using the first pattern that yields n_roots values."""
    for name, pat in _PATTERNS:
        matches = {int(m.group(1)): float(m.group(2)) for m in pat.finditer(log_text)}
        if len(matches) >= n_roots:
            energies = [matches[i] for i in sorted(matches)[:n_roots]]
            return name, energies
    return "none", []


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare Run A (Multistate=all) vs Run B (EFFE) MS-CASPT2 energies"
    )
    parser.add_argument("--log-a",    required=True, help="Run A log (Multistate=all)")
    parser.add_argument("--log-b",    required=True, help="Run B EFFE log")
    parser.add_argument("--n-roots",  type=int, default=8, help="Number of roots (default 8)")
    parser.add_argument("--threshold", type=float, default=1e-7,
                        help="Max allowed energy difference in Ha (default 1e-7)")
    args = parser.parse_args()

    log_a = open(args.log_a, errors="replace").read()
    log_b = open(args.log_b, errors="replace").read()

    pat_a, energies_a = parse_energies(log_a, args.n_roots)
    pat_b, energies_b = parse_energies(log_b, args.n_roots)

    print("=" * 65)
    print("  N2 EFFE Validation — Energy Comparison")
    print("=" * 65)
    print(f"  Run A ({pat_a}): {len(energies_a)} roots found")
    print(f"  Run B ({pat_b}): {len(energies_b)} roots found")
    print(f"  Threshold: {args.threshold:.0e} Ha")
    print()

    if not energies_a or not energies_b:
        print("ERROR: could not parse energies from one or both logs.")
        print("       Check that the runs completed (Happy landing) and")
        print("       adjust parsing patterns if needed.")
        sys.exit(2)

    if len(energies_a) != args.n_roots or len(energies_b) != args.n_roots:
        print(f"WARNING: expected {args.n_roots} roots, got A={len(energies_a)} B={len(energies_b)}")

    n = min(len(energies_a), len(energies_b))
    header = f"  {'Root':>4}  {'E_A (Ha)':>20}  {'E_B (Ha)':>20}  {'|ΔE| (Ha)':>14}  {'Status':>6}"
    print(header)
    print("  " + "-" * 63)

    all_pass = True
    for i in range(n):
        diff = abs(energies_a[i] - energies_b[i])
        status = "OK" if diff <= args.threshold else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"  {i+1:>4}  {energies_a[i]:>20.10f}  {energies_b[i]:>20.10f}  {diff:>14.2e}  {status:>6}")

    print()
    if all_pass:
        print("  RESULT: PASS — all energies match within threshold.")
        print("          EFFE coupling extraction is numerically equivalent")
        print("          to Multistate=all. Safe to proceed with root batching.")
    else:
        print("  RESULT: FAIL — one or more roots exceed threshold.")
        print("          Check coupling extraction pattern and H_eff assembly.")

    print("=" * 65)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
