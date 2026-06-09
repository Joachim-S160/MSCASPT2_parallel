#!/usr/bin/env python3
"""
Parse a Po2 SO-CASSI benchmark OpenMolcas log and write structured metrics.

Extracts:
  - Per-step wall time (SEWARD, SCF, each RASSCF, each CASPT2, RASSI)
  - Reference weights per root per spin block
  - SO-RASSI State 1 energy
  - Peak memory from /usr/bin/time -v output
  - Success/failure status

Output: benchmark_metrics.json + benchmark_summary.txt

Also stubs a 'future_adaptive' section listing roots with ref_weight < threshold
(configurable, default 0.70) to support future adaptive root-count reduction.
"""

import argparse
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_module_timings(log_text: str) -> dict[str, float]:
    """
    Return {section_label: wall_time_s} by scanning module start/stop markers.

    Actual OpenMolcas log format (confirmed from local log files):
      --- Start Module: rasscf at <datetime> ---
        ...
      --- Stop Module: rasscf at <datetime> /rc=_RC_ALL_IS_WELL_ ---
      --- Module rasscf spent N seconds ---
          Timing: Wall=12.23 User=11.68 System=0.42

    Crucially: Timing: appears AFTER --- Stop Module ---, so we track the
    most-recently-closed module when assigning timing values.

    Label assignment uses position-ordering — only Gateway prints the user
    title inside the module section; RASSCF/CASPT2 modules do NOT echo their
    section-specific Title= keywords in the module output. The cascade is
    always GATEWAY → SEWARD → RASSCF(Q) → CASPT2(Q) → RASSCF(T) →
    CASPT2(T) → RASSCF(S) → CASPT2(S) → RASSI.

    For modules that appear multiple times (rasscf, caspt2), the Nth
    occurrence is labelled using the ordered spin list [quintet, triplet, singlet].
    """
    # The expected ordered sequence of module occurrences in the log.
    # This maps (module_name, occurrence_index) → label.
    SPIN_ORDER = ["quintet", "triplet", "singlet"]
    SINGLE_MODULES = {"gateway", "seward", "scf", "rassi"}  # appear exactly once

    timings: dict[str, list[float]] = {}
    # occurrence counters for repeated modules
    module_counts: dict[str, int] = {"rasscf": 0, "caspt2": 0}
    current_label: str = ""
    last_label: str = ""  # label of the most recently closed module

    re_start  = re.compile(r'--- Start Module:\s+(\w+)\s+at', re.IGNORECASE)
    re_stop   = re.compile(r'--- Stop Module:\s+(\w+)\s+at', re.IGNORECASE)
    re_timing = re.compile(r'Timing:\s+Wall\s*=\s*([\d.]+)', re.IGNORECASE)

    for line in log_text.splitlines():
        # Module start — resolve label immediately from occurrence order
        m = re_start.search(line)
        if m:
            mod = m.group(1).lower()
            if mod in SINGLE_MODULES:
                current_label = mod
            elif mod in module_counts:
                idx = module_counts[mod]
                if idx < len(SPIN_ORDER):
                    spin = SPIN_ORDER[idx]
                    current_label = f"{mod}_{spin}"
                else:
                    current_label = f"{mod}_{idx}"
                module_counts[mod] += 1
            else:
                current_label = mod
            continue

        # Module stop — freeze current_label for use by the Timing: line
        m = re_stop.search(line)
        if m:
            last_label = current_label
            continue

        # Timing: appears after Stop Module — assign to last closed module
        m = re_timing.search(line)
        if m and last_label:
            timings.setdefault(last_label, []).append(float(m.group(1)))

    # Take the LAST timing per label (handles retries / repeated sections)
    return {label: vals[-1] for label, vals in timings.items() if vals}


def parse_reference_weights(log_text: str) -> dict[str, list[float]]:
    """
    Return {spin_label: [ref_weight_root1, ref_weight_root2, ...]} for each
    spin block.

    OpenMolcas Multistate=all output format:
      Reference weight:           0.97178    (standalone line inside CASPT2 section)

    Spin blocks are identified by the order of CASPT2 module sections:
      1st caspt2 = quintet, 2nd = triplet, 3rd = singlet.
    (Title= keywords are NOT echoed inside the CASPT2 module output, only by Gateway.)
    """
    SPIN_ORDER = ["quintet", "triplet", "singlet"]
    spin_weights: dict[str, list[float]] = {s: [] for s in SPIN_ORDER}

    caspt2_count = 0
    in_caspt2 = False

    re_start_caspt2 = re.compile(r'--- Start Module:\s+caspt2\s+at', re.IGNORECASE)
    re_stop_caspt2  = re.compile(r'--- Stop Module:\s+caspt2\s+at', re.IGNORECASE)
    re_refweight    = re.compile(r'Reference weight:\s+([\d.]+)', re.IGNORECASE)

    for line in log_text.splitlines():
        if re_start_caspt2.search(line):
            in_caspt2 = True
            continue
        if re_stop_caspt2.search(line):
            in_caspt2 = False
            caspt2_count += 1
            continue

        if in_caspt2 and caspt2_count < len(SPIN_ORDER):
            m = re_refweight.search(line)
            if m:
                spin_weights[SPIN_ORDER[caspt2_count]].append(float(m.group(1)))

    return {k: v for k, v in spin_weights.items() if v}


def parse_so_rassi_energies(log_text: str) -> list[float]:
    """Return list of SO-RASSI State energies (Ha), ordered by state number."""
    pattern = re.compile(
        r'SO-RASSI State\s+(\d+)\s+Total energy:\s+([-\d.]+)', re.IGNORECASE
    )
    results: dict[int, float] = {}
    for m in pattern.finditer(log_text):
        state = int(m.group(1))
        energy = float(m.group(2))
        results[state] = energy
    return [results[k] for k in sorted(results)]


def parse_peak_memory_mb(time_file_text: str) -> float | None:
    """Parse 'Maximum resident set size (kbytes): N' from /usr/bin/time -v output."""
    m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', time_file_text)
    if m:
        return round(int(m.group(1)) / 1024, 1)  # KB → MB
    return None


def parse_molcas_mem_peak(log_text: str) -> float | None:
    """Parse OpenMolcas 'Largest memory requirement' summary (MB or GB)."""
    # Pattern: "Largest memory requirement: XXX MB" or "XXX GB"
    m = re.search(r'Largest memory requirement[:\s]+([\d.]+)\s*(MB|GB)', log_text, re.IGNORECASE)
    if m:
        val = float(m.group(1))
        unit = m.group(2).upper()
        return round(val * 1024 if unit == "GB" else val, 1)
    return None


def flag_low_ref_weights(
    ref_weights: dict[str, list[float]],
    threshold: float = 0.70,
) -> dict[str, list[int]]:
    """Return {spin: [1-based root indices with ref_weight < threshold]}."""
    flagged: dict[str, list[int]] = {}
    for spin, weights in ref_weights.items():
        bad = [i + 1 for i, w in enumerate(weights) if w < threshold]
        if bad:
            flagged[spin] = bad
    return flagged


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Collect Po2 SO-CASSI benchmark metrics")
    parser.add_argument("--log",         required=True, help="OpenMolcas log file")
    parser.add_argument("--time-file",   default="",    help="/usr/bin/time -v output file")
    parser.add_argument("--total-wall",  type=float, default=0,
                        help="Total wall time in seconds (from PBS job timer)")
    parser.add_argument("--pbs-jobid",   default="",    help="PBS job ID")
    parser.add_argument("--n-elec",      type=int, default=0,  help="CAS active electrons")
    parser.add_argument("--n-orb",       type=int, default=0,  help="CAS active orbitals")
    parser.add_argument("--ciroot-q",    type=int, default=0,  help="Quintet CIROOT")
    parser.add_argument("--ciroot-t",    type=int, default=0,  help="Triplet CIROOT")
    parser.add_argument("--ciroot-s",    type=int, default=0,  help="Singlet CIROOT")
    parser.add_argument("--molcas-mem",  type=int, default=0,  help="MOLCAS_MEM (MB, requested)")
    parser.add_argument("--ref-weight-threshold", type=float, default=0.70,
                        help="Flag roots with ref_weight below this (default 0.70)")
    parser.add_argument("--output",  default="benchmark_metrics.json",
                        help="Output metrics JSON file")
    parser.add_argument("--summary", default="benchmark_summary.txt",
                        help="Output human-readable summary file")
    args = parser.parse_args()

    log_text = Path(args.log).read_text(errors="replace")
    time_text = Path(args.time_file).read_text(errors="replace") if args.time_file and Path(args.time_file).exists() else ""

    # Parse all metrics
    timings    = parse_module_timings(log_text)
    ref_weights = parse_reference_weights(log_text)
    so_energies = parse_so_rassi_energies(log_text)
    peak_rss_mb = parse_peak_memory_mb(time_text)
    molcas_peak_mb = parse_molcas_mem_peak(log_text)
    success     = bool(re.search(r'happy landing', log_text, re.IGNORECASE))
    flagged     = flag_low_ref_weights(ref_weights, args.ref_weight_threshold)

    # Compute derived totals
    step_order = [
        "gateway", "seward", "scf",
        "rasscf_quintet", "caspt2_quintet",
        "rasscf_triplet", "caspt2_triplet",
        "rasscf_singlet", "caspt2_singlet",
        "rassi",
    ]
    caspt2_total = sum(
        timings.get(k, 0) for k in ["caspt2_quintet", "caspt2_triplet", "caspt2_singlet"]
    )
    rasscf_total = sum(
        timings.get(k, 0) for k in ["rasscf_quintet", "rasscf_triplet", "rasscf_singlet"]
    )
    total_measured = sum(timings.get(k, 0) for k in step_order)

    # Build metrics dict
    metrics = {
        "molecule": "Po2",
        "geometry": "eq (R=2.797 Å)",
        "pbs_jobid": args.pbs_jobid,
        "cascade": "Q-T-S",
        "method": "Multistate-all",
        "roots_per_job": None,  # None = all roots in one job (this baseline)
        "cas": {
            "n_elec": args.n_elec,
            "n_orb":  args.n_orb,
        },
        "roots": {
            "quintet": args.ciroot_q,
            "triplet": args.ciroot_t,
            "singlet": args.ciroot_s,
            "total":   args.ciroot_q + args.ciroot_t + args.ciroot_s,
        },
        "success": success,
        "timing_s": {
            step: round(timings.get(step, 0), 2) for step in step_order
        } | {
            "caspt2_total":   round(caspt2_total, 2),
            "rasscf_total":   round(rasscf_total, 2),
            "total_measured": round(total_measured, 2),
            "total_wall":     round(args.total_wall, 2),
        },
        "memory_mb": {
            "molcas_mem_requested": args.molcas_mem,
            "peak_rss_mb":          peak_rss_mb,
            "molcas_log_peak_mb":   molcas_peak_mb,
        },
        "reference_weights": {
            spin: {
                "values": weights,
                "n_roots": len(weights),
                "min": round(min(weights), 5) if weights else None,
                "mean": round(sum(weights) / len(weights), 5) if weights else None,
                "n_below_threshold": len([w for w in weights if w < args.ref_weight_threshold]),
            }
            for spin, weights in ref_weights.items()
        },
        "adaptive_future": {
            "ref_weight_threshold": args.ref_weight_threshold,
            "roots_to_drop": flagged,
            "note": (
                "Roots listed under roots_to_drop have ref_weight < threshold. "
                "Future adaptive runs may reduce CIROOT to exclude these."
            ),
        },
        "so_rassi_energies_ha": so_energies[:20],  # first 20 states
        "so_rassi_state1_ha": so_energies[0] if so_energies else None,
    }

    Path(args.output).write_text(json.dumps(metrics, indent=2))
    print(f"Metrics written: {args.output}")

    # ---------------------------------------------------------------------------
    # Human-readable summary
    # ---------------------------------------------------------------------------
    lines = [
        "=" * 60,
        "  Po2 SO-CASSI Roots Benchmark — Summary",
        "=" * 60,
        f"  Geometry       : eq (R=2.797 Å, bundled po2_eq.xyz)",
        f"  CAS            : CAS({args.n_elec},{args.n_orb})",
        f"  Roots (Q/T/S)  : {args.ciroot_q}/{args.ciroot_t}/{args.ciroot_s}"
        f"  = {args.ciroot_q + args.ciroot_t + args.ciroot_s} total",
        f"  Method         : Multistate=all  (baseline)",
        f"  Status         : {'SUCCESS' if success else 'FAILED'}",
        "",
        "  --- Timing (wall seconds) ---",
    ]

    label_width = max(len(s) for s in step_order) + 4
    for step in step_order:
        t = timings.get(step, 0)
        bar = "#" * max(1, int(t / 30))  # 1 char ≈ 30 s
        lines.append(f"  {step:<{label_width}} {t:>8.1f} s  {bar}")

    lines += [
        "",
        f"  RASSCF total   : {rasscf_total:>8.1f} s",
        f"  CASPT2 total   : {caspt2_total:>8.1f} s",
        f"  Measured total : {total_measured:>8.1f} s",
        f"  PBS wall total : {args.total_wall:>8.1f} s",
        "",
        "  --- Memory ---",
        f"  MOLCAS_MEM (requested) : {args.molcas_mem} MB",
    ]
    if peak_rss_mb:
        lines.append(f"  Peak RSS (/usr/bin/time) : {peak_rss_mb} MB")
    if molcas_peak_mb:
        lines.append(f"  Peak (OpenMolcas log)    : {molcas_peak_mb} MB")

    lines += ["", "  --- Reference weights ---"]
    for spin in ("quintet", "triplet", "singlet"):
        if spin not in ref_weights:
            lines.append(f"  {spin:10s}: (not parsed)")
            continue
        info = metrics["reference_weights"][spin]
        bad  = len([w for w in ref_weights[spin] if w < args.ref_weight_threshold])
        lines.append(
            f"  {spin:10s}: n={info['n_roots']:3d}  min={info['min']:.5f}"
            f"  mean={info['mean']:.5f}"
            f"  n<{args.ref_weight_threshold:.2f}={bad}"
        )
    if flagged:
        lines += [
            "",
            f"  *** Roots below threshold ({args.ref_weight_threshold}):",
        ]
        for spin, roots in flagged.items():
            lines.append(f"      {spin}: roots {roots}")
        lines.append(
            "  *** Future adaptive run may reduce CIROOT to exclude these."
        )

    if so_energies:
        lines += [
            "",
            "  --- SO-RASSI energies ---",
            f"  State 1 : {so_energies[0]:.10f} Ha",
        ]
        if len(so_energies) >= 2:
            gap = (so_energies[1] - so_energies[0]) * 219474.63  # Ha → cm⁻¹
            lines.append(f"  State 2 : {so_energies[1]:.10f} Ha  (gap = {gap:.1f} cm⁻¹)")

    lines.append("=" * 60)
    summary_text = "\n".join(lines) + "\n"

    Path(args.summary).write_text(summary_text)
    print(summary_text)


if __name__ == "__main__":
    main()
