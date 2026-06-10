#!/usr/bin/env python3
"""
Po2 parallel SO-CASSI — H_eff assembler and KILL_INTRUDERS filter.

Reads all  worker_{SPIN}_k*_b*.log  files for a given spin block, extracts
the H_eff coupling rows and reference weights from the CASPT2(Only=N) output,
applies the KILL_INTRUDERS threshold (truncate at first failing root), and
writes:

  assemble_{SPIN}.input   OpenMolcas input for CASPT2(EFFE) → produces JobMix
  heff_{SPIN}.json        H_eff matrix + per-root ref weights + truncation metadata

The EFFE matrix is always the FULL n_roots × n_roots matrix — truncation is
applied later at the RASSI level via Nr of JobIphs = 3 N_Q N_T N_S.

Usage:
  cd SO_CASSI_test/po2_parallel_socassi/
  python3 collect_heff.py --spin triplet --n-roots 90 --kill-intruders 0.70
"""

from __future__ import annotations  # X | None and list[X] on Python 3.8/3.9

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# H_eff row extraction  (logic from n2_effe_validation/assemble_effe.py)
# ---------------------------------------------------------------------------

def extract_all_couplings(log_text: str, n_roots: int) -> np.ndarray:
    """
    Return H_eff[n_roots, n_roots] by scanning all ONLY=N coupling sections.

    Each CASPT2(Only=N) block writes one "Hamiltonian Effective Couplings"
    section whose N-th row contains the couplings to all other roots.
    """
    header = "hamiltonian effective couplings"
    text_lower = log_text.lower()

    positions: list[int] = []
    start = 0
    while True:
        pos = text_lower.find(header, start)
        if pos < 0:
            break
        positions.append(pos)
        start = pos + 1

    if len(positions) != n_roots:
        raise ValueError(
            f"Expected {n_roots} 'Hamiltonian Effective Couplings' sections, "
            f"found {len(positions)}. Check that all ONLY= blocks converged."
        )

    # Scientific-notation values only: e.g. -1.089E+02 (excludes timing lines)
    coupling_re = re.compile(r"<[^|]+\|\s+(-?[\d.]+[Ee][+-]?\d+)")

    H_eff = np.zeros((n_roots, n_roots))
    for i, pos in enumerate(positions):
        end = positions[i + 1] if i + 1 < len(positions) else len(log_text)
        section = log_text[pos:end]
        couplings = [float(m) for m in coupling_re.findall(section)]
        if len(couplings) != n_roots:
            raise ValueError(
                f"Root {i + 1}: expected {n_roots} coupling values, got {len(couplings)}."
            )
        H_eff[i, :] = couplings

    asymmetry = float(np.max(np.abs(H_eff - H_eff.T)))
    print(f"  H_eff assembled ({n_roots}x{n_roots}), max asymmetry: {asymmetry:.2e}")
    if asymmetry > 1e-6:
        print("  WARNING: H_eff asymmetry > 1e-6 — check coupling extraction.")

    return H_eff


# ---------------------------------------------------------------------------
# Reference weight extraction
# ---------------------------------------------------------------------------

# Confirmed format from actual Po2 CASPT2 log (one block per ONLY=N call):
#   FINAL CASPT2 RESULT:
#       ...
#       Reference weight:           0.95334
_REF_WEIGHT_RE = re.compile(r"Reference weight:\s+([\d.]+)")


def extract_ref_weights(log_text: str, n_roots: int) -> list[float]:
    """Return reference weights in root order (one per CASPT2 block)."""
    weights = [float(m) for m in _REF_WEIGHT_RE.findall(log_text)]
    if len(weights) != n_roots:
        raise ValueError(
            f"Expected {n_roots} reference weights, found {len(weights)}. "
            f"Check that all ONLY= blocks produced 'Reference weight:' lines."
        )
    return weights


# ---------------------------------------------------------------------------
# KILL_INTRUDERS: truncate at first failing root
# ---------------------------------------------------------------------------

def find_truncation_point(ref_weights: list[float], threshold: float) -> int:
    """
    Return 1-based index of first root below threshold; -1 if none fails.
    All roots from that index onward are considered intruder-contaminated.
    """
    for i, w in enumerate(ref_weights):
        if w < threshold:
            return i + 1  # 1-based
    return -1


# ---------------------------------------------------------------------------
# Log concatenation: gather all worker logs for a spin in root order
# ---------------------------------------------------------------------------

def load_spin_logs(spin: str, workdir: Path) -> tuple[str, list[tuple[int, int]]]:
    """
    Gather all worker_{spin}_k*_b*.log files, sort by ROOT_START inferred from
    batch index, and concatenate. Returns (combined_log_text, batch_ranges).

    Because ROOT_START is not embedded in the filename, we sort by batch index
    (b{N}) which corresponds to the submission order in submit_socassi.sh.
    """
    log_files = sorted(workdir.glob(f"worker_{spin}_k*_b*.log"))
    if not log_files:
        print(f"  ERROR: no worker_{spin}_k*_b*.log files found in {workdir}",
              file=sys.stderr)
        sys.exit(1)

    # Sort by batch index extracted from filename
    def _batch_idx(p: Path) -> int:
        m = re.search(r"_b(\d+)\.log$", p.name)
        return int(m.group(1)) if m else 0

    log_files = sorted(log_files, key=_batch_idx)

    combined = ""
    ranges: list[tuple[int, int]] = []
    for lf in log_files:
        text = lf.read_text(errors="replace")
        combined += text
        # Extract root range from the "Roots:" line in the header
        m_start = re.search(r"ROOT_START[^\d]*(\d+)", text, re.IGNORECASE)
        m_end   = re.search(r"ROOT_END[^\d]*(\d+)",   text, re.IGNORECASE)
        if m_start and m_end:
            ranges.append((int(m_start.group(1)), int(m_end.group(1))))
        else:
            # Fall back: count CASPT2 blocks to infer range
            n = len(re.findall(r"hamiltonian effective couplings", text, re.IGNORECASE))
            ranges.append((-1, n))

    print(f"  Loaded {len(log_files)} log file(s) for spin={spin}")
    return combined, ranges


# ---------------------------------------------------------------------------
# EFFE input writer
# ---------------------------------------------------------------------------

_SEWARD_FILES = (
    "RunFile", "OneInt", "SymInfo", "NqGrid",
    "ChVec1", "ChSel1", "ChRed", "ChDiag", "ChMap", "ChRst",
)

_SPIN_LABEL = {"quintet": "S5", "triplet": "S3", "singlet": "S1"}


def write_effe_input(
    spin: str,
    n_roots: int,
    H_eff: np.ndarray,
    output_path: Path,
) -> None:
    """
    Write assemble_{spin}.input — Seward injection + CASPT2(EFFE) block.

    GATEWAY+SEWARD are NOT included (run only once in run_setup.pbs).
    Seward files are injected from shared/ just like workers.
    """
    spin_label = _SPIN_LABEL[spin]

    injection = "\n".join(
        f">>COPY $CurrDir/shared/seward.{f}  $WorkDir/$Project.{f}"
        for f in _SEWARD_FILES
    )

    effe_rows = "\n".join(
        " ".join(f"{H_eff[i, j]:.14E}" for j in range(n_roots))
        for i in range(n_roots)
    )
    effe_block = f"{n_roots}\n{effe_rows}"

    content = f"""\
* Po2 parallel SO-CASSI — Phase 3 assembly: CASPT2(EFFE) for {spin}
* Seward files injected from shared/ — no GATEWAY/SEWARD/RASSCF re-run.
* H_eff matrix: {n_roots}x{n_roots} (full, KILL_INTRUDERS applied at RASSI level).
* JobMix is produced in MOLCAS_WORKDIR; run_assemble.pbs copies it to shared/.

{injection}
>>COPY $CurrDir/shared/rasscf_{spin}.JobIph  $WorkDir/$Project.JobIph

&CASPT2
  MAXITER=  300
  Frozen=   78
  Multistate= all
  Imaginary Shift= 0.25
  EFFE
{effe_block}
End of Input
"""
    output_path.write_text(content)
    print(f"  Written: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Assemble H_eff from parallel CASPT2 workers, apply KILL_INTRUDERS, "
                    "write CASPT2(EFFE) input for assembly job."
    )
    parser.add_argument("--spin",            required=True,
                        choices=["quintet", "triplet", "singlet"],
                        help="Spin block to process")
    parser.add_argument("--n-roots",         type=int, required=True,
                        help="Total roots in this spin block")
    parser.add_argument("--kill-intruders",  type=float, default=0.70,
                        help="Reference weight threshold (default: 0.70)")
    parser.add_argument("--workdir",         default=".",
                        help="Directory containing worker log files (default: .)")
    args = parser.parse_args()

    workdir   = Path(args.workdir).resolve()
    spin      = args.spin
    n_roots   = args.n_roots
    threshold = args.kill_intruders

    print(f"\n=== collect_heff.py: spin={spin}  n_roots={n_roots}  kill_intruders={threshold} ===")

    # --- Load and concatenate worker logs ---
    combined_log, batch_ranges = load_spin_logs(spin, workdir)

    # --- Extract H_eff matrix ---
    print(f"\n[1] Extracting H_eff rows...")
    H_eff = extract_all_couplings(combined_log, n_roots)

    # --- Extract reference weights ---
    print(f"\n[2] Extracting reference weights...")
    ref_weights = extract_ref_weights(combined_log, n_roots)

    min_w  = min(ref_weights)
    mean_w = sum(ref_weights) / len(ref_weights)
    print(f"  Min ref weight: {min_w:.5f}  Mean: {mean_w:.5f}")

    # --- Apply KILL_INTRUDERS ---
    trunc_at = find_truncation_point(ref_weights, threshold)
    if trunc_at == -1:
        n_surviving = n_roots
        print(f"  No roots below threshold {threshold} — all {n_roots} roots survive.")
    else:
        n_surviving = trunc_at - 1
        print(f"  First root below {threshold}: root {trunc_at}  "
              f"→ keeping roots 1–{n_surviving} ({n_surviving}/{n_roots})")

    # --- Write EFFE input (full H_eff, truncation handled by RASSI) ---
    print(f"\n[3] Writing EFFE input...")
    assemble_input = workdir / f"assemble_{spin}.input"
    write_effe_input(spin, n_roots, H_eff, assemble_input)

    # --- Write metadata JSON (consumed by run_rassi.pbs) ---
    metadata = {
        "spin": spin,
        "n_roots": n_roots,
        "n_surviving": n_surviving,
        "truncation_at": trunc_at,
        "kill_intruders_threshold": threshold,
        "ref_weights": ref_weights,
        "ref_weight_min": min_w,
        "ref_weight_mean": mean_w,
    }
    json_path = workdir / f"heff_{spin}.json"
    json_path.write_text(json.dumps(metadata, indent=2))
    print(f"  Written: {json_path}")

    # --- Summary ---
    print(f"\n=== Summary: {spin} ===")
    print(f"  H_eff:       {n_roots}x{n_roots}")
    print(f"  n_surviving: {n_surviving} (RASSI will use Nr of JobIphs ... {n_surviving} ...)")
    print(f"  Min Ref.wt.: {min_w:.5f}")
    if trunc_at != -1:
        dropped = [i + 1 for i, w in enumerate(ref_weights) if i + 1 >= trunc_at]
        print(f"  Dropped roots ({len(dropped)}): {dropped[:10]}{'...' if len(dropped) > 10 else ''}")
    print()


if __name__ == "__main__":
    main()
