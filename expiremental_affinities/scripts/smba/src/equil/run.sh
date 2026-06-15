#!/bin/bash --login
#SBATCH -C [nfh|nal|nif|nvf|neh]
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=4:0:0
#SBATCH --account=vermaaslab
#SBATCH --job-name=runDUMMY_NAME

cd $SLURM_SUBMIT_DIR

module use /mnt/home/vermaasj/modules
module load NAMD/3.0.1-gpu

NUM=$(ls system_run[0-9][0-9][0-9].coor 2>/dev/null | wc -l)
printf -v PRINTNUM "%03d" "$NUM"

# stop after 11 runs (000..010) bc we know only 11 max runs in this case
if [ "$NUM" -ge 11 ]; then
    echo "All runs already completed (NUM=$NUM). Exiting."
    exit 0
fi

echo "(output index $PRINTNUM)"
srun namd3 \
    +p12 +setcpuaffinity +idlepoll +devices 0 \
    system_run.namd \
    > system_run${PRINTNUM}.log \
    2> system_run${PRINTNUM}.err

