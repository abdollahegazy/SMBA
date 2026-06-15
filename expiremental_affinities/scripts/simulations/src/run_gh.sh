#!/bin/bash --login
#SBATCH --constraint=NOAUTO:grace
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --time=4:0:0
#SBATCH --account=vermaaslab
#SBATCH -J __JOB_NAME__

# Safety: bail if this sim's normal OR _gh sibling is already RUNNING (see run.sh).
BASE_NAME=${SLURM_JOB_NAME%_gh}
RUNNING=$(squeue --me -h -t RUNNING --format="%j %i" | grep -E "^${BASE_NAME}(_gh)? " | grep -v "$SLURM_JOB_ID" | wc -l)
if [ "$RUNNING" -gt 0 ]; then
    echo "Another job named $SLURM_JOB_NAME is already running, exiting to avoid race conditions"
    exit 0
fi

source /cvmfs/software.eessi.io/versions/2023.06/init/bash && module use /opt/modules/all
module --force purge

cd $SLURM_SUBMIT_DIR
NAMD=/mnt/home/vermaasj/namd/Linux-ARM64-g++/namd3

# --- benchmark 2000 steps to size the production run to the remaining walltime ---
cp system_run.namd _bench.namd
sed -i "s/^run .*/run 2000/" _bench.namd
sed -i "s/^dcdfreq.*/#&/" _bench.namd
sed -i "s/^restartfreq.*/#&/" _bench.namd
sed -i "s/^outputTiming.*/outputTiming 100/" _bench.namd
sed -i "s/^outputname.*/outputname _bench_out/" _bench.namd
srun $NAMD +p8 +setcpuaffinity +devices 0 _bench.namd > _bench.log 2>&1

NS_PER_DAY=$(grep "^TIMING" _bench.log | awk '{for(i=1;i<=NF;i++) if($i=="ns/days,") print $(i-1)}' | awk '{s+=$1; n++} END {print s/n}')
if [ -z "$NS_PER_DAY" ]; then
    echo "ERROR: Benchmark failed, dumping log:"; cat _bench.log
    rm -f _bench.namd _bench.log _bench_out.*; exit 1
fi
rm -f _bench.namd _bench.log _bench_out.*

# steps that fit in remaining walltime (80% safety, 30 min buffer); 2 fs timestep
TIME_LEFT_SEC=$(( SLURM_JOB_END_TIME - $(date +%s) - 1800 ))
NSTEPS=$(python3 -c "n=int(${NS_PER_DAY} * 0.80 / 86400.0 * ${TIME_LEFT_SEC} / 2e-6); print(n - n % 100)")
echo "Benchmark: ${NS_PER_DAY} ns/day, ${TIME_LEFT_SEC}s remaining -> ${NSTEPS} steps"

sed -i "s/^run .*/run $NSTEPS/" system_run.namd

NUM=$(ls system_run[0-9][0-9][0-9].dcd 2>/dev/null | wc -l)
printf -v PRINTNUM "%03d" "$NUM"
echo "(output index $PRINTNUM)"
srun $NAMD \
    +p8 +setcpuaffinity +devices 0 \
    system_run.namd \
    > system_run${PRINTNUM}.log \
    2> system_run${PRINTNUM}.err
