from pathlib import Path
import subprocess
import time
from tqdm import tqdm

INFERENCE_IN_PATH = Path("/mnt/scratch/hegazyab/project/dockingAF3/input")
INFERENCE_OUT_PATH = Path("/mnt/scratch/hegazyab/project/dockingAF3/output")
INFERENCE_OUT_PATH.mkdir(parents=True, exist_ok=True)

SPECIES = ["Arabidopsis"] #, "DouglasFir", "Eucalyptus", "Human"]
MAX_JOBS = 950
CHECK_INTERVAL = 20  # seconds

LOG_PATH = Path("/mnt/scratch/hegazyab/project/dockingAF3/logs")
LOG_PATH    .mkdir(parents=True, exist_ok=True)

def get_job_count(user="hegazyab"):
    result = subprocess.run(["squeue", "-u", user], capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    return max(len(lines) - 1, 0)  # subtract header

def submit_job(script_str):
    """Submit SLURM job using sbatch with the script passed as a string."""
    process = subprocess.Popen(["sbatch"], stdin=subprocess.PIPE, text=True)
    process.communicate(input=script_str)
    print("[SUBMITTED]")

# We'll keep a local counter of submitted jobs in this run.
local_job_count = 0

for s in tqdm(SPECIES, desc="Species"):
    print(f"\nWORKING ON SPECIES {s}")

    species_in = INFERENCE_IN_PATH / s
    species_out = INFERENCE_OUT_PATH / s
    species_out.mkdir(parents=True, exist_ok=True)

    species_logs = LOG_PATH / s
    species_logs.mkdir(parents=True,exist_ok=True)

    for protein_dir in species_in.iterdir():
        if not protein_dir.is_dir():
            continue

        output_dir_base = species_out / protein_dir.name
        output_dir_base.mkdir(parents=True, exist_ok=True)

        protein_logs = species_logs / protein_dir.name
        protein_logs.mkdir(parents=True, exist_ok=True)

        for ligand_json in protein_dir.iterdir():
            if ligand_json.suffix != ".json":
                continue

            ligand_out_dir = output_dir_base / ligand_json.stem 

            ligand_log_dir = protein_logs / ligand_json.stem

            expected_output = ligand_out_dir / f"{ligand_out_dir.stem}_model.cif"

            if expected_output.exists():
                 print(f"[SKIP] {expected_output} exists;")
                 continue


            job_name = f"{s}_{protein_dir.stem}_{ligand_json.stem}"

            docking_script = f"""#!/bin/bash --login
#SBATCH --account=vermaaslab
#SBATCH -C v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --job-name={job_name}
#SBATCH --mail-type=TIME_LIMIT
#SBATCH --mail-user=hegazyab@msu.edu
#SBATCH --output={ligand_log_dir}/slurm-%j.out
#SBATCH --error={ligand_log_dir}/slurm-%j.err

mkdir -p {ligand_log_dir}

singularity exec \\
    --nv \\
    --env XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" \\
    --bind {ligand_json}:/root/af_input/input.json \\
    --bind {output_dir_base}:/root/af_output \\
    --bind /mnt/home/hegazyab/af3:/root/model \\
    --bind /mnt/research/common-data/alphafold/database_3:/root/public_databases \\
    /mnt/home/hegazyab/af3/alphafold3.sif \\
    python /mnt/home/hegazyab/af3/run_alphafold.py \\
    --input_dir=/root/af_input \\
    --model_dir=/root/model \\
    --db_dir=/root/public_databases \\
    --output_dir=/root/af_output/ \\
    --flash_attention_implementation=xla \\
    --norun_data_pipeline
"""
            # If our local counter is below MAX_JOBS, we can submit without checking.
            if local_job_count >= MAX_JOBS:
                # Update local counter by checking squeue
                while get_job_count() >= MAX_JOBS:
                    print(f"[WAIT] Job count exceeds {MAX_JOBS}. Sleeping {CHECK_INTERVAL}s...")
                    time.sleep(CHECK_INTERVAL)
                local_job_count = get_job_count()

            submit_job(docking_script)
            local_job_count += 1
