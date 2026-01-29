#!/bin/bash
set -euo pipefail

BASE="../complexes"

find "$BASE" -name run.sh | while read -r run_script; do
    sim_dir=$(dirname "$run_script")

    echo ">>> Checking $sim_dir"

    ###########################################################
    # STEP 0 — detect MOL2 crash ONLY using .err
    ###########################################################
    crashed=0

    if find "$sim_dir" -maxdepth 1 -name "*.err" -size +0c | grep -q .; then
        crashed=1
        echo "    ❌ Crashed (.err non-empty)"
    fi

    ###########################################################
    # STEP 1 — Submit MOL2 if not running/queued and not crashed
    ###########################################################
    jobname=$(grep -E '^#SBATCH -J' "$sim_dir/run.sh" | awk '{print $3}')
    
    if [[ -z "$jobname" ]]; then
        echo "    ⚠️  Could not extract job name, skipping"
        echo
        continue
    fi

    ###########################################################
    # STEP 2 — Submit if not crashed and not already running
    ###########################################################
    if [[ $crashed -eq 1 ]]; then
        echo "    🔕 Not submitting (crashed)"
        echo
        continue
    fi

    # Check if already queued/running
    if squeue -u "$USER" -o "%.100j" --noheader | grep -Fxq "$jobname"; then
    	echo "    ⏩ Already queued/running: $jobname"
    	echo
    	continue
    fi

    # Submit the job
    echo "    🟢 Submitting: $jobname"
    ( cd "$sim_dir" && sbatch run.sh )

    echo
done
