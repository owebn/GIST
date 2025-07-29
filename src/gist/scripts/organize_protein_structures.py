import argparse
import pathlib
import shutil
import sys


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    for gene_dir in args.input_dir.iterdir():
        if not gene_dir.is_dir():
            continue

        for af_file in gene_dir.iterdir():
            if af_file.is_dir():
                continue

            if af_file.name == ".DS_Store":
                continue

        shutil.copy(af_file, args.input_dir.joinpath(f"{gene_dir.stem}.cif"))


if __name__ == "__main__":
    main(sys.argv[1:])
