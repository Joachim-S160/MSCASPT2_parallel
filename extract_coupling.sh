#!/bin/bash

# Check if a file or pattern is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <log_file_or_pattern>"
    exit 1
fi

# Input parameter
input=$1

# Extract effective couplings
grep -A6 -i 'Hamiltonian Effective Couplings' "$input" | grep '<' | awk '{print $NF}'
