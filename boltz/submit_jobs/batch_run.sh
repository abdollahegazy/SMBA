#!/bin/bash --login
#SBATCH --account=vermaaslab
#SBATCH -C [neh|nel|nal|nif|nvf]
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=3:50:00
#SBATCH --job-name=boltz-batch-%A_%a
#SBATCH --array=1-40
#SBATCH --output=logs/batch_%A_%a.out
#SBATCH --error=logs/batch_%A_%a.err


cd $SLURM_SUBMIT_DIR   

mkdir -p logs


# pick ligands for this array index
LIGAND_LIST=$SLURM_SUBMIT_DIR/batches/batch_${SLURM_ARRAY_TASK_ID}.txt


echo "Processing batch ${SLURM_ARRAY_TASK_ID} with $(wc -l < "$LIGAND_LIST") ligands"

while IFS= read -r ligand_dir; do
    if [[ -d "$ligand_dir" ]] && [[ -x "$ligand_dir/run.sh" ]]; then
        echo "Processing: $ligand_dir"
        (cd "$ligand_dir" && bash run.sh > "slurm_${SLURM_ARRAY_TASK_ID}_${RANDOM}.out" 2>&1)
    else
        echo "ERROR: Missing directory or run.sh in $ligand_dir" >&2
    fi
done < "$LIGAND_LIST"

echo "Completed batch ${SLURM_ARRAY_TASK_ID}"
