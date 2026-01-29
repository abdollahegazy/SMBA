#!/bin/bash --login
#SBATCH -C [nfh|neh|nal|nif|nvf]
#SBATCH --array=1-12%1
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:0:0
#SBATCH --account=vermaaslab
#SBATCH -J defaultjob

cd $SLURM_SUBMIT_DIR

module use /mnt/home/vermaasj/modules
module load NAMD/3.0.1-gpu

NUM=$(ls system_run[0-9][0-9][0-9].dcd 2>/dev/null | wc -l)
printf -v PRINTNUM "%03d" "$NUM"


echo "(output index $PRINTNUM)"
srun namd3 \
    +p8 +setcpuaffinity \
    +devices 0 \
    system_run.namd \
    > system_run${PRINTNUM}.out \
    2> system_run${PRINTNUM}.e
