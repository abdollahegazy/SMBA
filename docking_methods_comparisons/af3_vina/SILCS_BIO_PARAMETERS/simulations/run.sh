#!/bin/bash --login
#SBATCH -C [neh|nal|nif|nvf]
#SBATCH --array=1-9%1
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:0:0
#SBATCH --account=vermaaslab

cd $SLURM_SUBMIT_DIR

module use /mnt/home/vermaasj/modules
module load NAMD/3.0.1-gpu

NUM=$(ls system_run[0-9][0-9][0-9].coor 2>/dev/null | wc -l)
printf -v PRINTNUM "%03d" "$NUM"


echo "(output index $PRINTNUM)"
srun namd3 \
    +p4 +setcpuaffinity \
    +devices 0 \
    system_run.namd \
    > system_run${PRINTNUM}.out \
    2> system_run${PRINTNUM}.e


# ++ppn 1