#!/bin/bash
# Sync important files from HPC speedup_mscaspt2 folder to local repo.
# Usage: bash sync_from_hpc.sh

HPC_USER="vsc45694"
HPC_HOST="login.hpc.ugent.be"   # <-- change if needed
HPC_SRC="/dodrio/scratch/projects/starting_2025_097/speedup_mscaspt2/"
LOCAL_DST="$(dirname "$0")/"

rsync -av \
    --include="*/" \
    --include="*.py" --include="*.sh" --include="*.inp" \
    --include="*.xyz" --include="*.yml" --include="*.md" \
    --include="READ.me" --include="*.pbs" \
    --include="sorbaldehyde-04.RasOrb" --include="template_root.inp" \
    --exclude="mscaspt2_venv/" \
    --exclude="*.h5" --exclude="*.RICDLib" --exclude="*.molden*" \
    --exclude="*.GssOrb" --exclude="*.guessorb.*" --exclude="*.LprOrb" \
    --exclude="*.LoProp.*" --exclude="*.SpdOrb.*" --exclude="*.RasOrb.[0-9]*" \
    --exclude="*.log" --exclude="*.oldlog" --exclude="*.err" --exclude="*.status" \
    --exclude="pbs.out" --exclude="pbs.err" --exclude="workflow.out" --exclude="workflow.err" \
    --exclude="runinfo/" --exclude="test_root*/" --exclude="xmldump" \
    --exclude="*.save" --exclude="*.save.*" \
    "${HPC_USER}@${HPC_HOST}:${HPC_SRC}" "${LOCAL_DST}"
