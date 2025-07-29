import argparse
import csv
import pathlib
import re
import sys

from Bio import SeqIO
from pandas import read_csv


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("infile", type=pathlib.Path)
    parser.add_argument("reference", type=pathlib.Path)
    parser.add_argument("outfile", type=pathlib.Path)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    uniprot_format = re.compile("UniProtKB/Swiss-Prot")

    df = read_csv(args.infile)
    gene_set = set([gene.upper() for gene in df["Gene"]])
    record = SeqIO.read(args.reference, "genbank")

    outdata = list()

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = feature.qualifiers.get("gene", [""])[0]
        if gene.upper() not in gene_set:
            continue

        ID = gene.upper()
        sequence = feature.qualifiers["translation"][0]

        db_references = feature.qualifiers["db_xref"]
        uniprot = None

        for db_xref in db_references:
            match = re.match(uniprot_format, db_xref)
            if match is not None:
                uniprot = db_xref.split(":")[1]

        if uniprot is None:
            continue

        outdata.append({"id": ID,
                        "name": gene,
                        "sequence": sequence,
                        "uniprot": uniprot})
        fieldnames = ["id", "name", "sequence", "uniprot"]

    with open(args.outfile, "w") as filehandle:
        csv_writer = csv.DictWriter(filehandle, fieldnames=fieldnames)
        csv_writer.writeheader()

        for row in outdata:
            csv_writer.writerow(row)


if __name__ == "__main__":
    main(sys.argv[1:])
