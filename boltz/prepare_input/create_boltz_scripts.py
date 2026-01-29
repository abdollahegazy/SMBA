from pathlib import Path
from tqdm import tqdm

OUT_DIR = Path("../input") 

def create_run_scripts():
    """
    Iterates all ligand directories under OUT_DIR and creates a run.sh script
    next to each run.yaml.

    """

    yaml_files = list(OUT_DIR.rglob("run.yaml"))
    
    print(f"Found {len(yaml_files)} run.yaml files")

    for yaml_file in tqdm(yaml_files, desc="Creating run.sh scripts"):
        ligand_dir = yaml_file.parent
        run_sh = ligand_dir / "run.sh"

        script_content = f"""PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:64 boltz predict run.yaml --accelerator gpu --out_dir ./"""

        # if not run_sh.exists():
        run_sh.write_text(script_content)
        # else:
            # continue


if __name__ == "__main__":
    create_run_scripts()
