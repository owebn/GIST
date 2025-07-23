import argparse
import csv
import pathlib
import pickle
import shutil
import sys

from rdkit import Chem
import pymol

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


HEADER = ["pdbid", "-logKd/Ki"]
DIST = 4


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_dir", type=pathlib.Path,
                        required=True)
    parser.add_argument("-o", "--output_dir", type=pathlib.Path,
                        required=True)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True, parents=True)
    complexes_dir = output_dir.joinpath("complexes")
    complexes_dir.mkdir(exist_ok=True)
    metadata_file = output_dir.joinpath("metadata.csv")

    preprocess(args.input_dir, complexes_dir, metadata_file)


def preprocess(input_dir, complexes_dir, metadata_file):
    data_dicts = list()
    counter = 0
    for af_pred_cif in input_dir.iterdir():
        if not af_pred_cif.is_file():
            continue

        if not af_pred_cif.suffix == ".cif":
            continue

        if af_pred_cif.name == ".DS_Store":
            continue

        af_name = af_pred_cif.stem
        protein_name, ligand_name, _ = af_name.split("_")

        valid = generate_files(af_pred_cif, protein_name, ligand_name,
                               complexes_dir)
        if not valid:
            continue

        data_dicts.append({HEADER[0]: "".join([protein_name, ligand_name]),
                           HEADER[1]: counter})
        counter += 1

    write_metadata(data_dicts, metadata_file)


def write_metadata(data_dicts, metadata_file):
    header = HEADER
    with open(metadata_file, "w") as filehandle:
        csv_writer = csv.DictWriter(filehandle, fieldnames=header)

        csv_writer.writeheader()

        for row in data_dicts:
            csv_writer.writerow(row)


def generate_files(af_pred_cif, protein_name, ligand_name, complexes_dir,
                   distance=DIST):
    complex_name = "".join([protein_name, ligand_name])

    complex_dir = complexes_dir.joinpath(complex_name)
    complex_dir.mkdir(exist_ok=True)

    molecule_name = "_".join([protein_name, ligand_name])

    pocket_file = complex_dir.joinpath(f"{complex_name}_pocket.pdb")
    generate_pocket(af_pred_cif, pocket_file, molecule_name, ligand_name)

    ligand_file = complex_dir.joinpath(f"{complex_name}_ligand.pdb")
    generate_molecule(af_pred_cif, ligand_file, molecule_name, ligand_name)

    complex_path = complex_dir.joinpath(f"{complex_name}_{distance}A.rdkit")
    valid = generate_complex_file(pocket_file, ligand_file, complex_path)

    if not valid:
        shutil.rmtree(complex_dir)

    return valid


def generate_complex_file(pocket_file, ligand_file, complex_path):
    ligand = Chem.MolFromPDBFile(str(ligand_file.absolute()), removeHs=True)
    pocket = Chem.MolFromPDBFile(str(pocket_file.absolute()), removeHs=True)

    if pocket is None or ligand is None:
        if ligand is None:
            print(f"No molecule found for {ligand_file}")
        if pocket is None:
            print(f"No molecule found for {pocket_file}")

        return False

    complex_tuple = (ligand, pocket)
    with open(complex_path, "wb") as filehandle:
        pickle.dump(complex_tuple, filehandle)

    return True


def generate_pocket(prediction_path, pocket_path, molecule_name, ligand_name,
                    distance=5):
    ligand_chain_name = ligand_name.upper()
    pymol.cmd.load(str(prediction_path))

    pymol.cmd.remove("resn HOH")
    pymol.cmd.remove("hydrogens")

    chains = pymol.cmd.get_chains(molecule_name)
    assert ligand_chain_name in chains

    pymol.cmd.select(ligand_chain_name, f"chain {ligand_chain_name}")
    pymol.cmd.select("Pocket", f"byres {ligand_chain_name} around {distance}")
    pymol.cmd.save(str(pocket_path), "Pocket")
    pymol.cmd.delete("all")


def generate_molecule(prediction_path, chain_path, molecule_name, chain_name):
    chain_name = chain_name.upper()
    pymol.cmd.load(str(prediction_path))

    pymol.cmd.remove("resn HOH")
    pymol.cmd.remove("hydrogens")

    chains = pymol.cmd.get_chains(molecule_name)
    assert chain_name in chains

    pymol.cmd.select(chain_name, f"chain {chain_name}")
    pymol.cmd.save(str(chain_path), chain_name)
    pymol.cmd.delete("all")


if __name__ == "__main__":
    main(sys.argv[1:])
