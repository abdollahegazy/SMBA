#!/bin/bash
set -euo pipefail

BASE="../complexes/"

find "$BASE" -name "eq.sh" | while read -r eqsh; do
    dir=$(dirname "$eqsh")

    echo ">>> Checking $dir"

    (
        cd "$dir" || exit 1

        # ----------------------------------------------------
        # SKIP if system_eq001.coor already exists
        # ----------------------------------------------------
        if [[ -f system_eq001.coor ]]; then
            echo "⏭️  Skipping EQ — system_eq001.coor already exists"
            exit 0
        fi

        echo "    Submitting EQ job..."

        # Submit SLURM job
        set +e
        sbatch eq.sh
        rc=$?
        set -e

        # Report status
        if [[ $rc -ne 0 ]]; then
            echo "❌ [ERROR]: sbatch failed (exit $rc)"
        else
            echo "✅ EQ job submitted"
        fi
    )
    
    echo
done