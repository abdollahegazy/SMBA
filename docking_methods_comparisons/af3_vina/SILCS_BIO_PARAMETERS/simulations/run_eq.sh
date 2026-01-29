# #!/bin/ bash --login
# #SBATCH -C [neh|nel|nal|nif|nvf]
# #SBATCH --gres=gpu:1
# #SBATCH --gres-flags=enforce-binding
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --cpus-per-task=1
# #SBATCH --mem=16G
# #SBATCH --time=4:0:0
# #SBATCH --account=vermaaslab

# module use /mnt/home/vermaasj/modules
# module load NAMD/3.0.1-gpu


BASE_DIR="/mnt/scratch/hegazyab/dock_comp/simulations/complexes"
EQ_SCRIPT="eq.sh"


# # loop over every protein / ligand / {MD,AF3}_sim directory
find "$BASE_DIR" -type d \( -name "MD_sim" -o -name "AF3_sim" \) | while read -r simdir; do
    echo ">>> Running eq in $simdir"
    (
      cd "$simdir" && sbatch "$EQ_SCRIPT" 
    )
done

# echo "All jobs launched!"


    # # count how many previous .dcds to pick a new index
    # NUM=$(ls system_eq[0-9][0-9][0-9].dcd 2>/dev/null | wc -l)
    # #assign printnum to num in 3 digits
    # printf -v PRINTNUM "%03d" "$NUM"
    # # run the eq
    # srun namd3 ++ppn 1 +devices 0 system_eq.namd \
    #   > system_eq${PRINTNUM}.out \
    #   2> system_eq${PRINTNUM}.e
    # echo "  → launched with index $PRINTNUM in dir $simdir"