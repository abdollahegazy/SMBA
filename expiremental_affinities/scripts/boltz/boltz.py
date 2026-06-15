from utils import load_pairs, DATA_DIR

from varidock.stages import BoltzConfig, BoltzPredict
from varidock.types import BoltzInput, Ligand

#diffusionv2.py
# i had to move svd to cpu 
# i coudlnt figure out some linalg backend issue for my life lol

import os


# without this boltz launches a billion cpus and we reach the
# cpu walltime limit in like 4 min
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")


def make_boltz_input(protein_id: str, ligand_id: str, smiles: str) -> BoltzInput | None:
    # data_json = (
    #     DATA_DIR / protein_id / "af_output" / protein_id / f"{protein_id}_data.json"
    # )

    data_json = (
        DATA_DIR / protein_id /"MSA"/ protein_id / "af_output" / protein_id / f"{protein_id}_data.json"
    )

    output_dir = DATA_DIR / protein_id / "boltz" / ligand_id
    if not data_json.exists():
        print(f"Skipping {protein_id}: no AF3 data JSON")
        return None

    return BoltzInput(
        data_json_path=data_json,
        protein_chain_id="P",
        ligand=Ligand(name=ligand_id, smiles=smiles, af3_sequence_id="L"),
        output_dir=output_dir,
        name=f"{protein_id}_{ligand_id}",
    )


def run_one(args: tuple[BoltzInput, int]):
    inp, gpu_id = args
    predictor = BoltzPredict(BoltzConfig())
    # override env per process
    import os

    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    predictor.run(inp, write_only=False)


def is_boltz_done(protein_id: str, ligand_id: str) -> bool:
    name = f"{protein_id}_{ligand_id}"
    aff = (
        DATA_DIR
        / protein_id
        / "boltz"
        / ligand_id
        / "boltz_output"
        / f"boltz_results_{name}"
        / "predictions"
        / name
        / f"affinity_{name}.json"
    )
    return aff.exists()


def main():
    pairs = load_pairs()

    predictor = BoltzPredict(
        BoltzConfig()
    )

    for protein_id, ligand_id, smiles in pairs:
        if is_boltz_done(protein_id, ligand_id):
            print(f"Skipping {protein_id} + {ligand_id}: already done")
            continue

        inp = make_boltz_input(protein_id, ligand_id, smiles)
        if inp is None:
            print(f"Skipping {protein_id} + {ligand_id}: could not make BoltzInput")
            continue

        print(f"Running {protein_id} + {ligand_id}")
        try:
            predictor.run(inp, write_only=False)
        except Exception as e:
            print(f"Error running Boltz for {protein_id} + {ligand_id}: {e}")


if __name__ == "__main__":
    main()