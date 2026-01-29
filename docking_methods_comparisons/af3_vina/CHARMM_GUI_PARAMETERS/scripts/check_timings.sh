#!/bin/bash
set -euo pipefail

BASE="../complexes"

echo "=== Simulation Progress Report ==="
echo "Format: <nanoseconds>  <status>  <directory>"
echo

declare -A TIME_MAP
declare -A STATUS_MAP

# --- Helper: Replicate EXACT get_nanoseconds() logic ---
get_ns() {
    local dir="$1"

    # No .out → 2.020 ns
    if ! ls "$dir"/system_run*.out 1>/dev/null 2>&1; then
        echo "2.020"
        return
    fi

    # Latest out
    local latest_out=$(ls -1 "$dir"/system_run*.out | sort -V | tail -n 1)
    local last_step=$(grep -o "STEP [0-9]*" "$latest_out" | awk '{print $2}' | tail -n 1)

    # Fallback: second-to-last out
    if [[ -z "$last_step" ]]; then
        local second_last_out=$(ls -1 "$dir"/system_run*.out | sort -V | tail -n 2 | head -n 1)
        if [[ "$second_last_out" != "$latest_out" ]]; then
            last_step=$(grep -o "STEP [0-9]*" "$second_last_out" | awk '{print $2}' | tail -n 1)
        fi
    fi

    # Still no step → -1
    if [[ -z "$last_step" ]]; then
        echo "-1"
        return
    fi

    awk -v step="$last_step" 'BEGIN { printf "%.3f", step * 2e-6 }'
}

# --- Helper: detect crash ---
check_crashed() {
    local dir="$1"
    for efile in "$dir"/*.e; do
        if [[ -f "$efile" && -s "$efile" ]]; then
            if ! grep -q "CANCELLED AT" "$efile" 2>/dev/null; then
                return 0
            fi
        fi
    done
    return 1
}

# --- Collect ALL directories (mol2 + pdb) that contain timing ---
while IFS= read -r out; do
    dir=$(dirname "$out")
    TIME_MAP["$dir"]="$(get_ns "$dir")"

    if check_crashed "$dir"; then
        STATUS_MAP["$dir"]="CRASHED"
    fi
done < <(find "$BASE" -type f -name "system_run*.out")

# --- Now find all directories that SHOULD have timing but didn't appear ---
while IFS= read -r timing_dir; do
    if [[ -z "${TIME_MAP[$timing_dir]+x}" ]]; then
        TIME_MAP["$timing_dir"]="$(get_ns "$timing_dir")"

        if check_crashed "$timing_dir"; then
            STATUS_MAP["$timing_dir"]="CRASHED"
        fi
    fi
done < <(find "$BASE" -type d -regex ".*/\(pdb_params\|mol2_params\).*")

# --- Print sorted ---
for dir in "${!TIME_MAP[@]}"; do
    echo "${TIME_MAP[$dir]}  ${STATUS_MAP[$dir]:-}  $dir"
done | sort -n
