#!/bin/bash
# Configuration for the Po2 SO-CASSI roots benchmark.
# Sourced by run_baseline.pbs (and future run_effe_k.pbs variants).
#
# Cascade:  Quintet → Triplet → Singlet
# Start:    RHF SCF → RASSCF (no external RasOrb required)
# Method:   Multistate=all, HEFF   (Phase 1 baseline)
# Goal:     measure per-step wall time + peak RSS as reference before
#           benchmarking EFFE root-batching (Phase 2)

# ---------------------------------------------------------------------------
# Molecule — geometry bundled as po2_eq.xyz (R=2.797 Å ≈ Re)
# ---------------------------------------------------------------------------
MOL_LOWER="po2"
BASIS="ANO-RCC-VQZP"

# ---------------------------------------------------------------------------
# CAS — Po2 CAS(12,8), 168 total electrons
#   INACTIVE = (168 - 12) / 2 = 78
#   FROZEN   = 78  (same as INACTIVE for CASPT2)
# ---------------------------------------------------------------------------
CAS_N_ELEC=12
CAS_N_ORB=8
INACTIVE=78
FROZEN=78

# CIROOT values for the benchmark
CIROOT_Q=70
CIROOT_T=90
CIROOT_S=15

# ---------------------------------------------------------------------------
# CASPT2 settings
# ---------------------------------------------------------------------------
IMAGINARY_SHIFT="0.25"
LEVSHIFT="0.1"      # fixed — no retry loop in this benchmark

# ---------------------------------------------------------------------------
# Resource settings
# ---------------------------------------------------------------------------
MOLCAS_MEM=14000    # MB
OMP_NUM_THREADS=1

# ---------------------------------------------------------------------------
# HPC: path to setup_hortense.sh (needed to load modules + activate venv)
# ---------------------------------------------------------------------------
INSTALL_DIR="/dodrio/scratch/projects/starting_2025_097/autoCAS4HE_built/autoCAS4HE"
