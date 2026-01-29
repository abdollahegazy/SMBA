#!/usr/bin/env bash
set -euo pipefail

# assume this script lives alongside system_eq.namd
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# echo $SCRIPT_DIR
NAMD_CONF_1="$SCRIPT_DIR/system_eq.namd"
NAMD_CONF_RUN="$SCRIPT_DIR/system_run.namd"
EQ_SCRIPT="eq.sh"
RUN_SCRIPT="run.sh"

# base directory containing your complexes/
BASE_DIR="/mnt/scratch/hegazyab/dock_comp/simulations/complexes"

# # find each MD_sim and AF3_sim dir and copy the config in
find "$BASE_DIR" -type d \( -name "MD_sim" -o -name "AF3_sim" \) -print0 | \
while IFS= read -r -d '' simdir; do
  # cp "$NAMD_CONF_1" "$simdir/"
  # cp "$NAMD_CONF_RUN" "$simdir/"
  # cp "$EQ_SCRIPT" "$simdir/"
  cp "$RUN_SCRIPT" "$simdir/"
  # echo "Copied system_eq.namd and system_run.namd → $simdir"
done

echo "Done"


