#!/usr/bin/env python3
"""
Parse per-root H_eff couplings from input_B_roots log and write input_B_effe.input.

The B_roots log contains N_ROOTS sequential CASPT2(ONLY=N) blocks.  Each block
outputs a "Hamiltonian Effective Couplings" section with the N-th row of H_eff:

  Hamiltonian Effective Couplings
    < 1 |  -1.08965765213E+02     <- H_eff[N, 1]
    < 2 |  -1.09015273421E+02
    ...

Assembles the full N×N H_eff matrix and writes:

  &GATEWAY + &SEWARD + >>COPY rasscf_B.JobIph → &CASPT2 EFFE <matrix>
  + >>COPY JobMix.B + >>COPY JOB001 + &RASSI

Port of extract_couplings / create_combined_input from run_mscaspt2_Po2.py,
simplified for a single-spin validation with a sequential log.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np


def extract_all_couplings(log_text: str, n_roots: int) -> np.ndarray:
    """Return H_eff[n_roots, n_roots] by scanning all ONLY=N coupling sections."""
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

    # Scientific-notation values only: e.g. -1.089E+02
    # Fixed-decimal values (timing lines) are excluded by requiring [Ee].
    coupling_re = re.compile(r"<[^|]+\|\s+(-?[\d.]+[Ee][+-]?\d+)")

    H_eff = np.zeros((n_roots, n_roots))
    for i, pos in enumerate(positions):
        end = positions[i + 1] if i + 1 < len(positions) else len(log_text)
        section = log_text[pos:end]
        couplings = [float(m) for m in coupling_re.findall(section)]
        if len(couplings) != n_roots:
            raise ValueError(
                f"Root {i+1}: expected {n_roots} coupling values, got {len(couplings)}."
            )
        H_eff[i, :] = couplings

    # Symmetry check (informational)
    asymmetry = np.max(np.abs(H_eff - H_eff.T))
    print(f"H_eff assembled ({n_roots}x{n_roots}), max asymmetry: {asymmetry:.2e}")
    if asymmetry > 1e-6:
        print("WARNING: H_eff asymmetry > 1e-6 — check coupling extraction.")

    return H_eff


def write_effe_input(
    H_eff: np.ndarray,
    workdir: str,
    frozen: int = 2,
    imaginary_shift: float = 0.25,
) -> str:
    """Write input_B_effe.input to workdir and return its path."""
    n_roots = H_eff.shape[0]

    effe_block = "\n".join(
        [str(n_roots)]
        + [" ".join(f"{H_eff[i, j]:.14E}" for j in range(n_roots)) for i in range(n_roots)]
    )

    content = f"""\
* =============================================================================
* N2 EFFE validation — Run B phase 2 (EFFE combine)
* H_eff assembled from {n_roots} per-root ONLY= runs.
* JobIph injected from rasscf_B.JobIph (saved during B_roots run).
* =============================================================================

&GATEWAY
  Coord= n2.xyz
  Basis= cc-pVDZ
  Group= NoSym

&SEWARD
  Cholesky
End of Input

>>COPY $CurrDir/rasscf_B.JobIph $WorkDir/$Project.JobIph

&CASPT2
  MAXITER=  300
  Frozen=   {frozen}
  Multistate= all
  Imaginary Shift= {imaginary_shift}
  EFFE
{effe_block}

>>COPY $Project.JobMix $CurrDir/JobMix.B

>>COPY $CurrDir/JobMix.B JOB001
&RASSI
  Nr of JobIphs= 1 all
End of Input
"""
    out_path = Path(workdir) / "input_B_effe.input"
    out_path.write_text(content)
    print(f"Written: {out_path}")
    return str(out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble EFFE input from per-root B_roots log")
    parser.add_argument("--log",      required=True, help="Path to log_B_roots.log")
    parser.add_argument("--n-roots",  type=int, default=8, help="Number of roots (default 8)")
    parser.add_argument("--workdir",  required=True, help="Directory to write input_B_effe.input")
    parser.add_argument("--frozen",   type=int, default=2, help="Frozen orbitals for CASPT2 (default 2)")
    parser.add_argument("--imaginary-shift", type=float, default=0.25,
                        help="Imaginary shift for CASPT2 (default 0.25)")
    args = parser.parse_args()

    log_text = Path(args.log).read_text(errors="replace")
    H_eff = extract_all_couplings(log_text, args.n_roots)

    # Print H_eff diagonal for quick sanity check
    print("H_eff diagonal (root energies, Ha):")
    for i in range(args.n_roots):
        print(f"  Root {i+1:2d}: {H_eff[i, i]:.10f}")

    write_effe_input(
        H_eff,
        workdir=args.workdir,
        frozen=args.frozen,
        imaginary_shift=args.imaginary_shift,
    )


if __name__ == "__main__":
    main()
