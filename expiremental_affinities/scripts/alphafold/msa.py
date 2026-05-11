import os
import signal
from pathlib import Path
import pandas as pd
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict
from argparse import Namespace

from varidock.stages import AF3MSA, AF3InputBuilder
from varidock.runners.af3 import AF3Config
from varidock.types import ProteinSequence, AF3MSAOutput


def _worker_init():
    os.setpgrp()
    signal.signal(signal.SIGINT, signal.SIG_IGN)



def load_args() -> Namespace:
    parser = argparse.ArgumentParser(
        description="AF3 MSA generation for monomer proteins."
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Write inputs but don't execute AF3.",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run AF3 MSA jobs in parallel.",
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=7,
        help="Number of parallel workers to use if --parallel is set.",
    )
    return parser.parse_args()


def load_proteins() -> Dict[str, str]:
    data = pd.read_csv("../data/known_structures.csv")
    proteins = data['measurement_id'].map(lambda x: str(x).strip())
    seqs = data['sequences']
    return dict(zip(proteins, seqs))

def run_one(
    protein: ProteinSequence, builder: AF3InputBuilder, msa_stage: AF3MSA
) -> AF3MSAOutput:
    msa_input = builder.run(protein)

    expected = (
        msa_input.output_dir.resolve()
        / "af_output"
        / msa_input.protein_id
        / f"{msa_input.protein_id}_data.json"
    )


    if expected.exists():
        print(f"{protein.name} already done, skipping")
        return AF3MSAOutput(
            data_json_path=expected,
            protein_id=msa_input.protein_id,
            chain_id=msa_input.chain_id,
        )
    
    return msa_stage.run(msa_input)


if __name__ == "__main__":
    args = load_args()
    proteins = load_proteins()
    print(f"Loaded {len(proteins)} proteins")

    chain_id = "P"
    af3_cfg = AF3Config.from_config(script_args=("--norun_inference",))
    msa_stage = AF3MSA(af3_cfg, write_only=args.dry_run, overwrite_input=True)

    protein_inputs = [
        ProteinSequence(name=pid, sequence=seq) for pid, seq in proteins.items()
    ]

    if args.parallel:
        pool = ProcessPoolExecutor(max_workers=args.max_workers, initializer=_worker_init)
        futures = {}
        try:
            for p in protein_inputs:
                out_dir = Path(f"../data/predictions/{p.name}/MSA").resolve()
                out_dir.mkdir(parents=True, exist_ok=True)
                builder = AF3InputBuilder(output_dir=out_dir, chain_id=chain_id)
                futures[pool.submit(run_one, p, builder, msa_stage)] = p.name
            for f in as_completed(futures):
                print(f"{futures[f]} done")
                f.result()
        except KeyboardInterrupt:
            for fut in futures:
                fut.cancel()
            for pid in list(pool._processes):
                try:
                    os.kill(pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

        pool.shutdown(wait=False,cancel_futures=True)
    else:
        for p in protein_inputs:
            out_dir = Path(f"../data/predictions/{p.name}/MSA").resolve()
            out_dir.mkdir(parents=True, exist_ok=True)
            builder = AF3InputBuilder(output_dir=out_dir, chain_id=chain_id)
            msa_output = run_one(p, builder, msa_stage)
            print(f"{p.name} done → {msa_output.data_json_path}")
