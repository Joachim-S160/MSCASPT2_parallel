#!/usr/bin/env python3
"""
Po2 k-scaling benchmark results collector.

Reads all  worker_k{K}_b{IDX}.log  and  worker_k{K}_b{IDX}.time  files in
the current directory, aggregates per-k, and outputs:
  - Table to stdout
  - scaling_results.json
  - scaling_results.png  (if --plot given)

Usage:
  cd SO_CASSI_test/po2_scaling_test/
  python3 collect_scaling.py [--spin triplet] [--plot scaling_results.png]
"""

import argparse
import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Datetime parsing (inlined from po2_roots_benchmark/collect_metrics.py)
# ---------------------------------------------------------------------------
_DT_FORMATS = (
    "%a %b %d %I:%M:%S %p %Y",   # HPC: "Tue Jun  9 10:58:25 AM CEST 2026" after tz strip
    "%a %b  %d %I:%M:%S %p %Y",  # double-space day variant
    "%a %b %d %H:%M:%S %Y",      # local: "Tue Jun  9 15:28:50 2026"
    "%a %b  %d %H:%M:%S %Y",     # double-space day variant
)


def _parse_module_dt(dt_raw: str) -> datetime | None:
    dt_str = re.sub(r'\s+[A-Z]{2,5}\s*$', '', dt_raw.strip())
    for fmt in _DT_FORMATS:
        try:
            return datetime.strptime(dt_str, fmt)
        except ValueError:
            pass
    return None


# ---------------------------------------------------------------------------
# Per-batch parsers
# ---------------------------------------------------------------------------

def parse_batch_wall(log_text: str) -> float:
    """Sum of (Stop − Start) for all caspt2 module blocks in the log."""
    re_start = re.compile(r'--- Start Module:\s+caspt2\s+at\s+(.+?)\s*---', re.IGNORECASE)
    re_stop  = re.compile(r'--- Stop Module:\s+caspt2\s+at\s+(.+?)\s*(?:/rc=|---)', re.IGNORECASE)

    starts = [_parse_module_dt(m.group(1)) for m in re_start.finditer(log_text)]
    stops  = [_parse_module_dt(m.group(1)) for m in re_stop.finditer(log_text)]

    total = 0.0
    for s, e in zip(starts, stops):
        if s and e:
            total += (e - s).total_seconds()
    return total


def parse_rss_mb(time_text: str) -> float | None:
    m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', time_text)
    if m:
        return round(int(m.group(1)) / 1024, 1)
    return None


def parse_n_roots_completed(log_text: str) -> int:
    return len(re.findall(r'--- Start Module:\s+caspt2\s+at', log_text, re.IGNORECASE))


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def collect_results(spin: str, n_roots_total: int) -> dict[int, dict]:
    """
    Scan cwd for worker_k{K}_b{IDX}.log files and aggregate per k.
    Returns {k: {batches: [...], parallel_wall_s, max_rss_mb, ...}}
    """
    log_files = sorted(Path(".").glob("worker_k*_b*.log"))
    if not log_files:
        print("No worker_k*_b*.log files found.", file=sys.stderr)
        return {}

    raw: dict[int, list[dict]] = {}
    for lf in log_files:
        m = re.match(r"worker_k(\d+)_b(\d+)\.log", lf.name)
        if not m:
            continue
        k = int(m.group(1))
        b = int(m.group(2))
        tf = lf.with_suffix(".time")

        log_text  = lf.read_text(errors="replace")
        time_text = tf.read_text(errors="replace") if tf.exists() else ""

        wall  = parse_batch_wall(log_text)
        rss   = parse_rss_mb(time_text)
        n_completed = parse_n_roots_completed(log_text)
        success = bool(re.search(r'happy landing', log_text, re.IGNORECASE))

        raw.setdefault(k, []).append({
            "batch_idx": b,
            "wall_s": wall,
            "rss_mb": rss,
            "n_roots_completed": n_completed,
            "success": success,
        })

    results = {}
    for k in sorted(raw):
        batches = raw[k]
        walls = [b["wall_s"] for b in batches if b["success"]]
        rsss  = [b["rss_mb"] for b in batches if b["success"] and b["rss_mb"] is not None]

        if not walls:
            print(f"  k={k}: no successful batches.", file=sys.stderr)
            continue

        n_batches_full = math.ceil(n_roots_total / k)
        max_rss = max(rsss) if rsss else None

        results[k] = {
            "k": k,
            "n_measured": len(batches),
            "n_success": len(walls),
            "batch_walls_s": walls,
            "mean_wall_s": round(sum(walls) / len(walls), 1),
            "max_wall_s":  round(max(walls), 1),
            "proj_parallel_s": round(max(walls), 1),  # bottleneck batch = total HPC wall
            "rss_mb": rsss,
            "max_rss_mb": max_rss,
            "n_batches_full": n_batches_full,
            "total_node_mem_gb": round(n_batches_full * max_rss / 1024, 1) if max_rss else None,
        }

    return results


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
BASELINE_RASSCF_S = 667
BASELINE_RASSI_S  = 1866

def print_table(results: dict[int, dict]) -> None:
    print()
    print("=" * 88)
    print("  Po2 CAS(12,8)/ANO-RCC-VQZP — Triplet CASPT2 k-scaling  (90 roots, 30.9 s/root est.)")
    print("=" * 88)
    hdr = (f"  {'k':>4}  {'n_ok':>5}  {'mean_wall(s)':>12}  "
           f"{'proj_par(s)':>11}  {'max_RSS(MB)':>11}  "
           f"{'n_batches_full':>14}  {'total_mem(GB)':>13}")
    print(hdr)
    print("  " + "-" * 84)
    for k, r in sorted(results.items()):
        total_par = r["proj_parallel_s"] + BASELINE_RASSCF_S + BASELINE_RASSI_S
        rss_str   = f"{r['max_rss_mb']:.1f}" if r["max_rss_mb"] else "N/A"
        mem_str   = f"{r['total_node_mem_gb']:.1f}" if r["total_node_mem_gb"] else "N/A"
        print(f"  {k:>4}  {r['n_success']:>5}  {r['mean_wall_s']:>12.1f}  "
              f"{r['proj_parallel_s']:>11.1f}  {rss_str:>11}  "
              f"{r['n_batches_full']:>14}  {mem_str:>13}")
    print("=" * 88)
    print(f"  Baseline non-parallelizable overhead: RASSCF {BASELINE_RASSCF_S} s + RASSI {BASELINE_RASSI_S} s = {BASELINE_RASSCF_S + BASELINE_RASSI_S} s")
    print(f"  Full parallel wall = proj_par + {BASELINE_RASSCF_S + BASELINE_RASSI_S} s")
    print()


def make_plot(results: dict[int, dict], output_path: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plot.", file=sys.stderr)
        return

    ks   = sorted(results)
    wall = [results[k]["proj_parallel_s"] for k in ks]
    rss  = [results[k]["max_rss_mb"] or 0 for k in ks]
    tmem = [results[k]["total_node_mem_gb"] or 0 for k in ks]

    fig, ax1 = plt.subplots(figsize=(7, 5))

    color_wall = "#2166ac"
    color_rss  = "#d6604d"
    color_mem  = "#f4a582"

    ax1.plot(ks, wall, "o-", color=color_wall, linewidth=2, markersize=7, label="proj. parallel wall (s)")
    ax1.axhline(BASELINE_RASSCF_S + BASELINE_RASSI_S, color="grey", linestyle="--",
                linewidth=1.2, label=f"non-parallelizable floor ({BASELINE_RASSCF_S + BASELINE_RASSI_S} s)")
    ax1.set_xlabel("k  (roots per worker)", fontsize=12)
    ax1.set_ylabel("Projected parallel wall (s)", color=color_wall, fontsize=11)
    ax1.tick_params(axis="y", labelcolor=color_wall)
    ax1.set_xticks(ks)

    ax2 = ax1.twinx()
    ax2.plot(ks, rss, "s--", color=color_rss, linewidth=2, markersize=7, label="peak RSS / worker (MB)")
    ax2.plot(ks, [t * 1024 for t in tmem], "^:", color=color_mem, linewidth=1.5,
             markersize=6, label="total node mem (MB, full sweep)")
    ax2.set_ylabel("Memory (MB)", color=color_rss, fontsize=11)
    ax2.tick_params(axis="y", labelcolor=color_rss)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=9)

    plt.title("Po2 triplet CASPT2 — k-scaling (90 roots, ANO-RCC-VQZP)", fontsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Plot saved: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Collect Po2 k-scaling benchmark results")
    parser.add_argument("--spin", default="triplet", help="Spin block (for metadata; default: triplet)")
    parser.add_argument("--n-roots", type=int, default=90, help="Total roots in spin block (default: 90)")
    parser.add_argument("--plot", default="", metavar="PNG", help="Save matplotlib graph to this path")
    parser.add_argument("--output", default="scaling_results.json", help="JSON output file")
    args = parser.parse_args()

    results = collect_results(args.spin, args.n_roots)
    if not results:
        print("No results found.", file=sys.stderr)
        sys.exit(1)

    print_table(results)

    # Write JSON
    payload = {
        "molecule": "Po2",
        "spin": args.spin,
        "n_roots_total": args.n_roots,
        "baseline_rasscf_s": BASELINE_RASSCF_S,
        "baseline_rassi_s": BASELINE_RASSI_S,
        "results": {str(k): v for k, v in results.items()},
    }
    Path(args.output).write_text(json.dumps(payload, indent=2))
    print(f"JSON written: {args.output}")

    if args.plot:
        make_plot(results, args.plot)


if __name__ == "__main__":
    main()
