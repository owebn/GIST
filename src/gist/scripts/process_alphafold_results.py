import argparse
import csv
import pathlib
import shutil
import sys


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("input_dir", type=pathlib.Path)
    parser.add_argument("output_dir", type=pathlib.Path)
    parser.add_argument("--metadata_file", type=pathlib.Path)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    args.output_dir.mkdir(exist_ok=True, parents=True)
    results_map = process_alphafold_results(args.input_dir, args.output_dir)

    if args.metadata_file is not None:
        outdicts = []
        for job_id, files in results_map.items():
            cif_file, summary_file = files

            outdicts.append({"id": job_id,
                             "cif_path": str(cif_file),
                             "summary_path": str(summary_file)})

        fieldnames = ["id", "cif_path", "summary_path"]
        with open(args.metadata_file, "w") as filehandle:
            csv_writer = csv.DictWriter(filehandle, fieldnames=fieldnames)

            csv_writer.writeheader()
            for row in outdicts:
                csv_writer.writerow(row)


def process_alphafold_results(input_dir, output_dir):
    results_map = dict()
    for alphafold_dir in input_dir.iterdir():
        if not alphafold_dir.is_dir():
            continue

        new_cif_file, new_summary_file = process_alphafold_result_dir(
                                                    alphafold_dir, output_dir)

        results_map[alphafold_dir.stem] = (new_cif_file, new_summary_file)

        shutil.rmtree(alphafold_dir)  # Consider removing?

    return results_map


def process_alphafold_result_dir(af_dir, output_dir):
    cif_name = f"{af_dir.stem}_model.cif"
    cif_file = af_dir.joinpath(cif_name)
    new_cif_file = output_dir.joinpath(cif_name)

    shutil.copy(cif_file, new_cif_file)

    summary_name = f"{af_dir.stem}_summary_confidences.json"
    summary_file = af_dir.joinpath(summary_name)
    new_summary_file = output_dir.joinpath(summary_name)

    shutil.copy(summary_file, new_summary_file)

    return new_cif_file, new_summary_file


if __name__ == "__main__":
    main(sys.argv[1:])
