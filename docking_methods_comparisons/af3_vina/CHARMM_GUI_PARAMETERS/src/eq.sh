#!/bin/bash --login
#SBATCH -C [neh|nal|nif|nvf]
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --account=vermaaslab
#SBATCH -J defaultjob

module use /mnt/home/vermaasj/modules
module load NAMD/3.0.1-gpu

cd "$SLURM_SUBMIT_DIR"  

runnum=$SLURM_ARRAY_TASK_ID

# pick the next free index for outputs
NUM=$(ls system_eq[0-9][0-9][0-9].dcd 2>/dev/null | wc -l)
printf -v IDX "%03d" "$NUM"

echo "runnum=$runnum  (output index $IDX)"
srun namd3 \
     +p4 +setcpuaffinity \
     +devices 0 \
     system_eq.namd \
  > system_eq${IDX}.out \
  2> system_eq${IDX}.err
