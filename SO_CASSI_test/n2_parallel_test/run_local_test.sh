#!/bin/bash
# =============================================================================
# N2 parallel CASPT2 local test
#
# Validates that GATEWAY + SEWARD + RASSCF can be run ONCE and the resulting
# files injected into worker jobs that run &CASPT2 Only=N truly in parallel.
#
# Workflow:
#   Phase 1:  GATEWAY + SEWARD + RASSCF  →  ./shared/ (11 files)
#   Seq ref:  inject 11 files + CASPT2(Only=1) + CASPT2(Only=2)  (one call)
#   Parallel: inject 11 files + CASPT2(Only=1)  (worker 1, background)
#             inject 11 files + CASPT2(Only=2)  (worker 2, background)
#   Compare:  row 1 and row 2 of H_eff from seq vs parallel, threshold 1e-7 Ha
#
# Key invariant: all pymolcas calls are cd'd to BENCH_DIR so $CurrDir = BENCH_DIR
# in every >>COPY directive.
#
# Usage:
#   cd SO_CASSI_test/n2_parallel_test/
#   bash run_local_test.sh
# =============================================================================
set -euo pipefail

BENCH_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export MOLCAS=/home/joaschee/OpenMolcas/build
PYMOLCAS=/home/joaschee/OpenMolcas/build/pymolcas
source /home/joaschee/molcas_env/bin/activate

export OMP_NUM_THREADS=1
export MOLCAS_MEM=4000
export MOLCAS_PRINT=verbose

cd "${BENCH_DIR}"
mkdir -p shared

echo "========================================"
echo "  N2 Parallel CASPT2 Local Test"
echo "========================================"
echo "Start: $(date)"
echo ""

# ---------------------------------------------------------------------------
# Phase 1: GATEWAY + SEWARD + RASSCF (run once, saves 11 files to ./shared/)
# ---------------------------------------------------------------------------
echo "=== PHASE 1: GATEWAY + SEWARD + RASSCF ==="
rm -rf mw_setup && mkdir -p mw_setup
MOLCAS_WORKDIR="${BENCH_DIR}/mw_setup" "${PYMOLCAS}" \
    "${BENCH_DIR}/phase1_setup.input" > "${BENCH_DIR}/phase1.log" 2>&1
if grep -qi "happy landing" "${BENCH_DIR}/phase1.log"; then
    echo "Phase 1: OK"
    echo "  Shared files written:"
    ls -lh "${BENCH_DIR}/shared/"
else
    echo "Phase 1: FAILED — check phase1.log"
    exit 1
fi

# ---------------------------------------------------------------------------
# Sequential reference: inject 11 files + CASPT2(Only=1) + CASPT2(Only=2)
# ---------------------------------------------------------------------------
echo ""
echo "=== SEQUENTIAL REFERENCE (roots 1 and 2) ==="
rm -rf mw_seq && mkdir -p mw_seq
MOLCAS_WORKDIR="${BENCH_DIR}/mw_seq" "${PYMOLCAS}" \
    "${BENCH_DIR}/phase2_seq.input" > "${BENCH_DIR}/phase2_seq.log" 2>&1
if grep -qi "happy landing" "${BENCH_DIR}/phase2_seq.log"; then
    echo "Seq ref: OK"
else
    echo "Seq ref: FAILED — check phase2_seq.log"
    exit 1
fi

# ---------------------------------------------------------------------------
# Parallel workers: CASPT2(Only=1) and CASPT2(Only=2) at the same time
# ---------------------------------------------------------------------------
echo ""
echo "=== PARALLEL WORKERS (roots 1 and 2 simultaneously) ==="
rm -rf mw_w1 mw_w2 && mkdir -p mw_w1 mw_w2
MOLCAS_WORKDIR="${BENCH_DIR}/mw_w1" "${PYMOLCAS}" \
    "${BENCH_DIR}/phase2_only1.input" > "${BENCH_DIR}/phase2_w1.log" 2>&1 &
PID1=$!
MOLCAS_WORKDIR="${BENCH_DIR}/mw_w2" "${PYMOLCAS}" \
    "${BENCH_DIR}/phase2_only2.input" > "${BENCH_DIR}/phase2_w2.log" 2>&1 &
PID2=$!

wait $PID1; RC1=$?
wait $PID2; RC2=$?

if grep -qi "happy landing" "${BENCH_DIR}/phase2_w1.log"; then
    echo "Worker 1: OK"
else
    echo "Worker 1: FAILED (exit ${RC1}) — check phase2_w1.log"
fi
if grep -qi "happy landing" "${BENCH_DIR}/phase2_w2.log"; then
    echo "Worker 2: OK"
else
    echo "Worker 2: FAILED (exit ${RC2}) — check phase2_w2.log"
fi

# ---------------------------------------------------------------------------
# Compare H_eff coupling rows: seq vs parallel
# ---------------------------------------------------------------------------
echo ""
echo "=== COMPARE COUPLING ROWS ==="
python3 - <<'PYEOF'
import re, sys

def get_coupling_rows(log_path):
    """Return list of coupling row lists, one per 'Hamiltonian Effective Couplings' section."""
    try:
        text = open(log_path, errors="replace").read()
    except FileNotFoundError:
        print(f"  ERROR: {log_path} not found")
        return []
    header = "hamiltonian effective couplings"
    low = text.lower()
    positions = []
    start = 0
    while True:
        pos = low.find(header, start)
        if pos < 0:
            break
        positions.append(pos)
        start = pos + 1
    rows = []
    for i, pos in enumerate(positions):
        end = positions[i + 1] if i + 1 < len(positions) else len(text)
        section = text[pos:end]
        vals = re.findall(r"<[^|]+\|\s+(-?[\d.]+[Ee][+-]?\d+)", section)
        rows.append([float(v) for v in vals])
    return rows

seq_rows = get_coupling_rows("phase2_seq.log")
w1_rows  = get_coupling_rows("phase2_w1.log")
w2_rows  = get_coupling_rows("phase2_w2.log")

print(f"  Sequential: {len(seq_rows)} coupling sections found")
print(f"  Worker 1:   {len(w1_rows)} coupling sections found")
print(f"  Worker 2:   {len(w2_rows)} coupling sections found")
print()

ok = True
for label, a_rows, a_idx, b_rows, b_idx in [
    ("Row 1  seq[0] vs w1[0]", seq_rows, 0, w1_rows, 0),
    ("Row 2  seq[1] vs w2[0]", seq_rows, 1, w2_rows, 0),
]:
    if a_idx >= len(a_rows) or b_idx >= len(b_rows):
        print(f"  {label}: PARSE ERROR")
        ok = False
        continue
    a, b = a_rows[a_idx], b_rows[b_idx]
    if len(a) != len(b):
        print(f"  {label}: LENGTH MISMATCH ({len(a)} vs {len(b)})")
        ok = False
        continue
    diffs = [abs(x - y) for x, y in zip(a, b)]
    maxd = max(diffs)
    status = "PASS" if maxd <= 1e-7 else "FAIL"
    if status == "FAIL":
        ok = False
    print(f"  {label}: max |Δ| = {maxd:.2e}  {status}")

print()
if ok:
    print("  RESULT: PASS — parallel workers reproduce sequential coupling rows.")
    print("          GATEWAY/SEWARD/RASSCF file injection is valid.")
    print("          True parallelism of CASPT2(ONLY=N) is confirmed.")
else:
    print("  RESULT: FAIL — check coupling extraction and logs above.")

sys.exit(0 if ok else 1)
PYEOF

echo ""
echo "========================================"
echo "  Done: $(date)"
echo "========================================"
