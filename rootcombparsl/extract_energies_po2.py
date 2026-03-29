#!/usr/bin/env python3
"""
Post-processing: collect SO-RASSI energies from all 25 Po2 geometry runs.

Reads:  CURVE_DIR/geom_NNN/final_rassi/final_rassi.log  (SO-RASSI energies)
        CURVE_DIR/geom_NNN/config_geom_NNN.yml           (xyz path for bond distance)

Writes: CURVE_DIR/mscaspt2_sorassi_10states.dat
        Format matches autoCAS4HE/tests/molcas/SOCASSI/VQZP_25/sorassi_10states.dat
"""

import re
import sys
from pathlib import Path

CURVE_DIR = Path("/dodrio/scratch/projects/2025_060/Joachim/po2_curve")
N_GEOM = 25
N_STATES = 10


def bond_distance_from_xyz(xyz_file: str) -> float:
    with open(xyz_file) as f:
        lines = f.readlines()
    x1 = float(lines[2].split()[1])
    x2 = float(lines[3].split()[1])
    return abs(x2 - x1)


def xyz_path_from_config(config_file: Path) -> str:
    with open(config_file) as f:
        for line in f:
            if 'xyz_file:' in line:
                return line.split('"')[1]
    raise ValueError(f"xyz_file not found in {config_file}")


def sorassi_energies(log_file: Path, n_states: int) -> list:
    content = log_file.read_text()
    matches = re.findall(
        r'SO-RASSI State\s+\d+\s+Total energy:\s+([-\d.]+)', content
    )
    return [float(e) for e in matches[:n_states]]


rows = []
missing = []

for i in range(N_GEOM):
    idx = f"{i:03d}"
    geom_dir = CURVE_DIR / f"geom_{idx}"
    rassi_log = geom_dir / "final_rassi" / "final_rassi.log"
    config = geom_dir / f"config_geom_{idx}.yml"

    if not rassi_log.exists():
        missing.append(idx)
        continue

    xyz = xyz_path_from_config(config)
    dist = bond_distance_from_xyz(xyz)
    energies = sorassi_energies(rassi_log, N_STATES)

    if len(energies) < N_STATES:
        print(f"WARNING: geom_{idx}: only {len(energies)}/{N_STATES} states found",
              file=sys.stderr)
        energies += [float('nan')] * (N_STATES - len(energies))

    rows.append((dist, energies))

if missing:
    print(f"Incomplete/missing geometries ({len(missing)}): {', '.join(missing)}",
          file=sys.stderr)

rows.sort(key=lambda r: r[0])

out_file = CURVE_DIR / "mscaspt2_sorassi_10states.dat"
with open(out_file, 'w') as f:
    header = "# Distance(A)  " + "  ".join(f"S{j+1}" for j in range(N_STATES))
    f.write(header + "\n")
    for dist, energies in rows:
        cols = "  ".join(f"{e:.6f}" for e in energies)
        f.write(f"{dist:.4f}  {cols}\n")

print(f"Written {len(rows)}/{N_GEOM} geometries to {out_file}")
