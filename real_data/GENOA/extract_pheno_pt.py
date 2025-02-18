#! /usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--file", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--keep", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

dd = pd.read_csv(f"{args.file}", sep="\t")
subj = pd.read_csv(f"{args.keep}", sep="\t", header=None)
dd = dd[["IID", args.gene]]
dd = dd[dd.IID.isin(subj[0])]

dd.to_csv(f"{args.out}.sushie.bed", sep="\t", header=False, index=False)
