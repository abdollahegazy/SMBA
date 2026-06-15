from argparse import ArgumentParser, Namespace

from varidock.runners.af3 import AF3Config
from varidock.stages import AF3Inference, AF3MSAMerger, AF3MSAMergerConfig
from varidock.types import AF3MergedInput, AF3MSAOutput, Ligand

from utils import DATA_DIR, load_pairs

PROTEIN_CHAIN_ID = "P"
LIGAND_CHAIN_ID = "L"


def load_args() -> Namespace:
    parser = ArgumentParser(
        description="Run AF3 protein+ligand inference for each pair in known_structures.csv."
    )
    parser.add_argument("--write-only", action="store_true")
    parser.add_argument("--overwrite-inp", action="store_true")
    return parser.parse_args()


def get_msa_output(protein_id: str) -> AF3MSAOutput:
    data_json = (
        DATA_DIR
        / protein_id
        / "MSA"
        / protein_id
        / "af_output"
        / protein_id
        / f"{protein_id}_data.json"
    )
    return AF3MSAOutput(
        data_json_path=data_json,
        protein_id=protein_id,
        chain_id=PROTEIN_CHAIN_ID,
    )


def is_inference_done(protein_id: str, ligand_id: str) -> bool:
    cif = (
        DATA_DIR
        / protein_id
        / "inference_af3_lig"
        / ligand_id
        / "af_output"
        / ligand_id
        / f"{ligand_id}_model.cif"
    )
    return cif.exists()


def run_one(
    af3_inference: AF3Inference,
    protein_id: str,
    ligand_id: str,
    smiles: str,
) -> None:
    msa = get_msa_output(protein_id)
    if not msa.data_json_path.exists():
        print(f"Skipping {protein_id}+{ligand_id}: MSA data JSON not found")
        return

    if is_inference_done(protein_id, ligand_id):
        print(f"{protein_id}+{ligand_id} already done, skipping")
        return

    inference_dir = DATA_DIR / protein_id / "inference_af3_lig" / ligand_id

    merger = AF3MSAMerger(
        AF3MSAMergerConfig(output_dir=DATA_DIR / protein_id / "inference_af3_lig")
    )
    merged = merger.run(
        msa_outputs=[msa],
        ligands=[
            Ligand(
                name=ligand_id,
                smiles=smiles,
                af3_sequence_id=LIGAND_CHAIN_ID,
            )
        ],
        name=ligand_id,
    )

    inference_input = AF3MergedInput(
        json_path=merged.json_path,
        name=ligand_id,
        output_dir=inference_dir,
    )

    print(f"Running AF3 inference for {protein_id}+{ligand_id}")
    af3_inference.run(inference_input)


def main():
    args = load_args()
    pairs = load_pairs()

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

    for protein_id, ligand_id, smiles in pairs:
        run_one(af3_inference, protein_id, ligand_id, smiles)


if __name__ == "__main__":
    main()
