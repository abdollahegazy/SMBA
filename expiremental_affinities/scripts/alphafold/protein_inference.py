import shutil
from argparse import ArgumentParser, Namespace
from pathlib import Path

from varidock.runners import AF3Config
from varidock.stages import AF3Inference
from varidock.types import AF3MergedInput


from utils import DATA_DIR, load_proteins



def load_args() -> Namespace:
    parser = ArgumentParser(description="Run AF3 inference on generated MSAs.")
    parser.add_argument("--write-only", action="store_true")
    parser.add_argument("--overwrite-inp", action="store_true")
    return parser.parse_args()


def run_one(af3_inference: AF3Inference, protein_id: str) -> None:
    if is_inference_done(protein_id):
        print(f"{protein_id} already done, skipping")
        return
    
    data_json = get_data_json_path(protein_id)

    if not data_json.exists():
        print(f"Skipping {protein_id}: no data JSON")
        return

    inference_dir = DATA_DIR / protein_id / "inference_af3"
    input_dir = inference_dir / "af_input"
    input_dir.mkdir(parents=True, exist_ok=True)

    dest_json = input_dir / f"{protein_id}.json"
    if not dest_json.exists():
        shutil.copy2(data_json, dest_json)

    inp = AF3MergedInput(
        json_path=dest_json,
        name=protein_id,
        output_dir=inference_dir,
    )

    print(f"Running inference for {protein_id}")
    af3_inference.run(inp)

def get_data_json_path(protein_id:str) -> Path:
    return  (
        DATA_DIR
        / protein_id
        / "MSA"
        / protein_id
        / "af_output"
        / protein_id
        / f"{protein_id}_data.json"
    )

def is_inference_done(protein_id: str) -> bool:
    cif = (
        DATA_DIR
        / protein_id
        / "inference_af3"
        / "af_output"
        / protein_id
        / f"{protein_id}_model.cif"
    )
    return cif.exists()



def main():
    args = load_args()
    proteins = sorted(load_proteins(DATA_DIR))


    af3_cfg = AF3Config.from_config(
        singularity_args=("--nv",),
        script_args=("--gpu_device=0",),
    )

    af3_inference = AF3Inference(
        af3_cfg,
        jax_cache_dir="/mnt/scratch/hegazyab/.cache",
        write_only=args.write_only,
        overwrite_input=args.overwrite_inp,
    )

    for protein_id in proteins:
        run_one(af3_inference, protein_id)


if __name__ == "__main__":
    main()
