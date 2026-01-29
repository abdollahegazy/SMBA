#!/bin/bash
set -euo pipefail

BASE="../complexes"

find "$BASE" -name system_eq.namd | while read -r eqfile; do
    dir=$(dirname "$eqfile")

    echo ">>> Checking $dir"

    (
        cd "$dir" || exit 1

        # ----------------------------------------------------
        # SKIP if system_eq000.coor already exists
        # ----------------------------------------------------
        if [[ -f system_eq000.coor ]]; then
            echo "⏭️  Skipping EQ — system_eq000.coor already exists"
            exit 0
        fi

        echo "Running EQ..."

        # Run NAMD EQ step
        set +e
        namd3 +p8 +setcpuaffinity +devices 1 system_eq.namd \
            > system_eq000.out \
            2> system_eq000.err
        rc=$?
        set -e

        # Report status
        if [[ $rc -ne 0 ]]; then
            echo "❌ [ERROR]: NAMD EQ failed (exit $rc)"
        else
            echo "✅ EQ Success"
        fi
    )
    
    echo
done
