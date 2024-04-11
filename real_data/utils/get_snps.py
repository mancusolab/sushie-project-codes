#! /usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--file", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

df = pd.read_csv(args.file, sep="\t")
df.snp.to_csv(args.out, sep="\t", index=False)

