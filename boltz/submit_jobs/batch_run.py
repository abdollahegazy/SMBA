from pathlib import Path
import subprocess
import time
from tqdm import tqdm
import hashlib

# Paths
INFERENCE_IN_PATH = Path("/mnt/scratch/hegazyab/boltz/input")
LOG_PATH = Path("/mnt/scratch/hegazyab/boltz/logs")

LOG_PATH.mkdir(parents=True, exist_ok=True)

# Config
MAX_JOBS = 980         # under cluster limit of 1000
CHECK_INTERVAL = 30    # seconds
USER = "hegazyab"

def get_job_count(user=USER):
    """Return number of jobs currently submitted/running for user."""
    result = subprocess.run(["squeue", "-u", user], capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    return max(len(lines) - 1, 0)  # subtract header

def submit_job(script_str):
    """Submit SLURM job using sbatch with the script passed as a string."""
    process = subprocess.run(["sbatch", "--parsable"], input=script_str,
                             text=True, capture_output=True)
    if process.returncode == 0:
        job_id = process.stdout.strip()
        print(f"[SUBMITTED] {job_id}")
    else:
        print(f"[ERROR submitting] {process.stderr}")

# Track jobs submitted this run
local_job_count = 0

# Iterate over all run.yaml files
yaml_files = list(INFERENCE_IN_PATH.rglob("run.yaml"))
print(f"Found {len(yaml_files)} run.yaml files")

for yaml_file in tqdm(yaml_files, desc="Submitting jobs"):
    ligand_dir = yaml_file.parent
    protein_dir = ligand_dir.parent
    species_dir = protein_dir.parent

    # expected output (customize to match what your run produces)
    expected_output = INFERENCE_IN_PATH / species_dir.name / protein_dir.name / ligand_dir.name / "boltz_results_run"
    if expected_output.exists():
        print(f"[SKIP] {expected_output} exists")
        continue

    # logs
    ligand_log_dir = LOG_PATH / species_dir.name / protein_dir.name / ligand_dir.name
    ligand_log_dir.mkdir(parents=True, exist_ok=True)

    # unique + short job name
    job_hash = hashlib.md5(str(ligand_dir).encode()).hexdigest()[:6]
    job_name = f"{species_dir.name}_{protein_dir.name}_{ligand_dir.name}_{job_hash}"

    run_sh = ligand_dir / "run.sh"
    if not run_sh.exists():
        print(f"[WARN] No run.sh found in {ligand_dir}, skipping")
        continue

    # Slurm job script
    job_script = f"""#!/bin/bash --login
#SBATCH --account=vermaaslab
#SBATCH -C [neh|nel|nal|nif|nvf]
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --job-name={job_name}
#SBATCH --output={ligand_log_dir}/slurm-%j.out
#SBATCH --error={ligand_log_dir}/slurm-%j.err

cd {ligand_dir}
bash run.sh &> log.run
"""

    # Throttle if too many jobs
    if local_job_count >= MAX_JOBS:
        while get_job_count() >= MAX_JOBS:
            print(f"[WAIT] Job count exceeds {MAX_JOBS}. Sleeping {CHECK_INTERVAL}s...")
            time.sleep(CHECK_INTERVAL)
        local_job_count = get_job_count()

    submit_job(job_script)
    local_job_count += 1
