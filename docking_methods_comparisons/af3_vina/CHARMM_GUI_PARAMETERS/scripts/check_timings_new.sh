#!/bin/bash
set -euo pipefail

BASE="../complexes"
MAX_JOBS=64

echo "=== Simulation Progress Report ==="
echo "Format: <nanoseconds_pdb> <status_pdb> | <nanoseconds_mol2> <status_mol2> | <path>"
echo

# --- Helper: Replicate EXACT get_nanoseconds() logic ---
get_ns() {
    local dir="$1"

    if ! ls "$dir"/system_run*.out 1>/dev/null 2>&1; then
        echo "2.020"
        return
    fi

    local latest_out=$(ls -1 "$dir"/system_run*.out | sort -V | tail -n 1)
    local last_step=$(grep -o "STEP [0-9]*" "$latest_out" | awk '{print $2}' | tail -n 1)

    if [[ -z "$last_step" ]]; then
        local second_last_out=$(ls -1 "$dir"/system_run*.out | sort -V | tail -n 2 | head -n 1)
        if [[ "$second_last_out" != "$latest_out" ]]; then
            last_step=$(grep -o "STEP [0-9]*" "$second_last_out" | awk '{print $2}' | tail -n 1)
        fi
    fi

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

# --- Helper: check if directory matches the target pattern ---
is_valid_dir() {
    local dir="$1"
    if [[ "$dir" =~ ../complexes/(pdb|mol2)_params/(boltz|vina)_lower/[^/]+/[^/]+/[^/]+/(boltz|vina)$ ]]; then
        return 0
    fi
    return 1
}

# --- Process a single directory ---
process_dir() {
    local dir="$1"
    local tmpfile="$2"
    
    if ! is_valid_dir "$dir"; then
        return
    fi
    
    local ns=$(get_ns "$dir")
    local status=""
    
    if check_crashed "$dir"; then
        status="CRASHED"
    fi
    
    echo "${dir}|${ns}|${status}" >> "$tmpfile"
}

# Create temp file for results
TMPFILE=$(mktemp)
trap "rm -f $TMPFILE" EXIT

# Get all directories
mapfile -t DIRS < <(find "$BASE" -type d -regex ".*/\(pdb_params\|mol2_params\)/\(boltz\|vina\)_lower/.*/.*/.*/\(boltz\|vina\)")

# Process in parallel with job control
for dir in "${DIRS[@]}"; do
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.1
    done
    
    process_dir "$dir" "$TMPFILE" &
done

wait

# Process results and pair pdb/mol2
declare -A PDB_DATA
declare -A MOL2_DATA

while IFS='|' read -r dir ns status; do
    # Extract the normalized path (everything after pdb_params or mol2_params)
    if [[ "$dir" =~ pdb_params/(.*) ]]; then
        normalized="${BASH_REMATCH[1]}"
        PDB_DATA["$normalized"]="${ns}  ${status}"
    elif [[ "$dir" =~ mol2_params/(.*) ]]; then
        normalized="${BASH_REMATCH[1]}"
        MOL2_DATA["$normalized"]="${ns}  ${status}"
    fi
done < "$TMPFILE"

# Combine and print
{
    for path in "${!PDB_DATA[@]}"; do
        pdb_info="${PDB_DATA[$path]}"
        mol2_info="${MOL2_DATA[$path]:-N/A}"
        echo "${pdb_info}|${mol2_info}|${path}"
    done

    # Also print mol2-only entries
    for path in "${!MOL2_DATA[@]}"; do
        if [[ -z "${PDB_DATA[$path]+x}" ]]; then
            mol2_info="${MOL2_DATA[$path]}"
            echo "N/A|${mol2_info}|${path}"
        fi
    done
} | sort -t'|' -k1 -n | awk -F'|' '{printf "%-20s | %-20s | %s\n", $1, $2, $3}'