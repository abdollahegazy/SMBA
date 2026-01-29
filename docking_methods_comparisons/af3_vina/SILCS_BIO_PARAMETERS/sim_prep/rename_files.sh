cd /path/to/complexes

find . -type d \( -name MD_sim -o -name AF3_sim \) | while read simdir; do
  for ext in psf pdb xsc; do
    # look for the file ending in _system.$ext
    src=$(ls "$simdir"/*_system.$ext 2>/dev/null)
    if [[ -f $src ]]; then
      mv -v "$src" "$simdir/system.$ext"
    fi
  done
done