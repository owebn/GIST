import argparse
import csv
import json
import pathlib
import sys


SEEDS = [1, 2]
DIALECT = "alphafold3"
VERSION = 1


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_proteins", type=pathlib.Path)
    parser.add_argument("input_ligands", type=pathlib.Path)
    parser.add_argument("outdir", type=pathlib.Path)
    parser.add_argument("--cif_dir", type=pathlib.Path, required=False)
    parser.add_argument("--drug_dir", action="store_true")

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    proteins = read_in_file(args.input_proteins)
    ligands = read_in_file(args.input_ligands)

    jobs = create_jobs(proteins, ligands, cif_dir=args.cif_dir)

    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    for job_id, job in jobs.items():
        if args.drug_dir:
            drug_id = job_id.split("_")[1]
            drug_dir = outdir.joinpath(drug_id)
            drug_dir.mkdir(exist_ok=True, parents=True)
                  
            write_json(drug_dir.joinpath(f"{job_id}.json"), job)

            continue

        write_json(outdir.joinpath(f"{job_id}.json"), job)


def create_jobs(proteins, ligands, cif_dir=None):
    jobs = dict()
    for protein in proteins:
        for ligand in ligands:
            protein_id = protein["id"]
            sequence = protein["sequence"]
            ligand_id = ligand["id"]
            smiles = ligand["smiles"]
            cif_path = None
            if cif_dir is not None:
                cif_path = str(cif_dir.joinpath(f"{protein_id}.cif"))

            job_id = f"{protein_id}_{ligand_id}"
            jobs[job_id] = create_job(job_id, protein_id, sequence,
                                      ligand_id, smiles, cif_path)

    return jobs


def create_job(job_id, protein_id, sequence, ligand_id, smiles, cif_path):
    indicies = list(range(len(sequence)))

    protein_dict = {"id": protein_id,
                    "sequence": sequence}
    ligand_dict = {"id": ligand_id,
                   "smiles": smiles.replace("\\", "\\")}

    if cif_path is not None:
        protein_dict["templates"] = [{"mmcifPath": cif_path,
                                      "queryIndices": indicies,
                                      "templateIndices": indicies}]
        protein_dict["pairedMsa"] = ""
        protein_dict["unpairedMsa"] = ""

    job_dict = {"name": job_id,
                "modelSeeds": SEEDS,
                "sequences": [
                    {"protein": protein_dict},
                    {"ligand": ligand_dict}],
                "dialect": DIALECT,
                "version": VERSION}

    return job_dict


def read_in_file(input_list):
    with open(input_list, "r") as filehandle:
        csv_reader = csv.DictReader(filehandle)

        rows = [row for row in csv_reader]

    return rows


def write_json(outfile, job_dict):
    with open(outfile, "w") as filehandle:
        filehandle.write(json.dumps(job_dict, indent=4, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv[1:])
