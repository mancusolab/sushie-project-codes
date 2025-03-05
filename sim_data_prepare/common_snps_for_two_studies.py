#! /usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--file1", type=str)
parser.add_argument("--file2", type=str)
parser.add_argument("--file3", type=str)
parser.add_argument("--file4", type=str)
parser.add_argument("--file5", type=str)

args = parser.parse_args()

df1 = pd.read_csv(args.file1, sep="\t", header=None)
df2 = pd.read_csv(args.file2, sep="\t", header=None)
df3 = pd.read_csv(args.file3, sep="\t", header=None)
df4 = pd.read_csv(args.file4, sep="\t", header=None)

n_row = df1[[1]].merge(df2[[1]], how="inner", on=1).merge(df3[[1]], how="inner", on=1).merge(df4[[1]], how="inner", on=1)
print(n_row.shape[0])
n_row.to_csv(args.file5, sep="\t", header=None, index=False)
