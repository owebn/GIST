import argparse
import pathlib
import sys
import urllib

from pandas import read_csv
from Bio.PDB.alphafold_db import get_structural_models_for


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("infile", type=pathlib.Path)
    parser.add_argument("outdir", type=pathlib.Path)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    df = read_csv(args.infile)

    for protein_id, uniprot_id in zip(df["id"], df["uniprot"]):
        protein_dir = outdir.joinpath(protein_id)
        protein_dir.mkdir(exist_ok=True)

        try:
            _ = [result for result in get_structural_models_for(
                                        uniprot_id, directory=protein_dir)]
        except urllib.error.HTTPError:
            continue


if __name__ == "__main__":
    main(sys.argv[1:])
