#! /usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--file", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--subject", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

dd = pd.read_csv(f"{args.file}", sep="\t")
subj = pd.read_csv(f"{args.subject}", sep="\t", header=None)
dd = dd[["geno_id", args.gene]]
dd = dd[dd.geno_id.isin(subj[0])]

dd.to_csv(f"{args.out}", sep="\t", header=False, index=False)
