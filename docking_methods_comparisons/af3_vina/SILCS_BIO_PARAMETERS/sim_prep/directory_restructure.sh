#!/usr/bin/env bash
cd ../complexes
set -euo pipefail

# adjust if your script lives somewhere else
BASE_DIR="./"

find "$BASE_DIR" -mindepth 3 -maxdepth 3 -type d | while read -r ligand_dir; do
  echo "→ Organizing $ligand_dir"

  # make target subdirs
  mkdir -p "$ligand_dir/MD_sim" "$ligand_dir/AF3_sim"

  # move MD_complex_sim_system.* into MD_sim
  for ext in psf pdb xsc; do
    src="$ligand_dir/MD_complex_sim_system.$ext"
    if [[ -e $src ]]; then
      mv "$src" "$ligand_dir/MD_sim/"
    fi
  done

  # move AF3_complex_sim_system.* into AF3_sim
  for ext in psf pdb xsc; do
    src="$ligand_dir/AF3_complex_sim_system.$ext"
    if [[ -e $src ]]; then
      mv "$src" "$ligand_dir/AF3_sim/"
    fi
  done

  # now remove every other regular file at this level
  find "$ligand_dir" -maxdepth 1 -type f -exec rm -f {} +

done

echo "Done!"
