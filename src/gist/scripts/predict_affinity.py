import argparse
import pathlib
import sys

import torch
import numpy as np
import pandas as pd

from dataset_GIGN import GraphDataset, PLIDataLoader
from utils import load_model_dict
from GIGN import GIGN

GRAPH_TYPE = "Graph_GIGN"
BATCH_SIZE = 128
DIST = 5


def parse_args(unparsed_args):
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_dir", type=pathlib.Path,
                        required=True)
    parser.add_argument("-o", "--output", type=pathlib.Path,
                        required=True)
    parser.add_argument("-m", "--model", type=pathlib.Path,
                        required=True)
    parser.add_argument("-np", "--cores", type=int, default=1)

    parser.add_argument("-g", "--graph_type", type=str, default=GRAPH_TYPE)
    parser.add_argument("-s", "--batch_size", type=str, default=BATCH_SIZE)

    args = parser.parse_args()
    return args


def main(unparsed_args):
    args = parse_args(unparsed_args)

    device = torch.device("cuda:0")
    initialize_dataloader(args.input_dir, args.cores,
                          graph_type=args.graph_type,
                          batch_size=args.batch_size,
                          create=True)

    dataloader = initialize_dataloader(args.input_dir, args.cores,
                                       graph_type=args.graph_type,
                                       batch_size=args.batch_size,
                                       create=False)
    model = load_model(device, args.model)

    pred, label = predict(model, dataloader, device)

    write_output(pred, label, args.output)


def write_output(pred, label, output_file):
    with open(output_file, "w") as filehandle:
        filehandle.write("pred,label\n")

        for x, y in zip(pred, label):
            filehandle.write(f"{x},{y}\n")


def initialize_dataloader(input_dir, workers,
                          graph_type=GRAPH_TYPE, batch_size=BATCH_SIZE,
                          dis_threshold=DIST, create=False):
    df_path = input_dir.joinpath("metadata.csv")
    dataset_path = input_dir.joinpath("complexes")

    df = pd.read_csv(df_path)
    dataset = GraphDataset(dataset_path, df, graph_type=graph_type,
                           dis_threshold=DIST,
                           create=create)
    dataloader = PLIDataLoader(dataset, batch_size=batch_size,
                               shuffle=create, num_workers=workers)

    return dataloader


def load_model(device, model_path):
    model = GIGN(35, 256).to(device)
    load_model_dict(model, (model_path))
    model = model.cuda()

    return model


def predict(model, dataloader, device):
    model.eval()

    pred_list = []
    label_list = []
    for data in dataloader:
        data = data.to(device)
        with torch.no_grad():
            pred = model(data)
            label = data.y

            pred_list.append(pred.detach().cpu().numpy())
            label_list.append(label.detach().cpu().numpy())

    pred = np.concatenate(pred_list, axis=0)
    label = np.concatenate(label_list, axis=0)

    data_df = dataloader.dataset.data_df
    data_map = dict()
    for cid, pka in zip(data_df["pdbid"], data_df["-logKd/Ki"]):
        data_map[pka] = cid

    label = [data_map[y] for y in label]

    return pred, label


if __name__ == "__main__":
    main(sys.argv[1:])
