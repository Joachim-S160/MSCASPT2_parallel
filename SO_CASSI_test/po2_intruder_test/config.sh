#!/bin/bash
# Po2 intruder-state validation test — shared configuration.
# Sourced by all PBS scripts and submit_intruder_test.sh.
#
# Root counts are DELIBERATELY ABOVE the Weyl-safe cap (62) for quintet
# to force intruder contamination at dissociation geometries. This tests
# the global KILL_INTRUDERS consensus mechanism.
#
#   Weyl max:  Q=70  T=378  S=336
#   Safe cap:  Q=62  T=62   S=62
#   This test: Q=50  T=62   S=25   ← Q in danger zone: roots ~44-50 expected to fail
#
# KILL_INTRUDERS=0.80 (stricter than production 0.70; easier to observe filter working)

# Root counts
CIROOT_Q=50
CIROOT_T=62
CIROOT_S=25

# Batching and filtering
K=5
KILL_INTRUDERS=0.80

# CAS parameters — Po2 CAS(12,8), ANO-RCC-VQZP
N_ELEC=12
N_ORB=8
INACTIVE=78
FROZEN=78
BASIS="ANO-RCC-VQZP"
IMAGINARY_SHIFT=0.25

# 8 representative geometries (R in Å: 1.95, 2.55, 2.80, 2.95, 3.20, 3.60, 5.00, 10.00)
N_GEOM=8
GEOM_LIST=(po2_000 po2_003 po2_006 po2_009 po2_012 po2_014 po2_018 po2_031)
# R values matching above (for refwt_all.dat distance column)
GEOM_R=(1.95 2.55 2.80 2.95 3.20 3.60 5.00 10.00)

INSTALL_DIR="/dodrio/scratch/projects/starting_2025_097/autoCAS4HE_built/autoCAS4HE"
