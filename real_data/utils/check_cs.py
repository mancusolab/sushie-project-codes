#! /usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--cs-path", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--cs-out", type=str)

args = parser.parse_args()

df_normal = pd.read_csv(f"{args.cs_path}/{args.gene}.normal.sushie.cs.tsv", sep="\t")
df_indep = pd.read_csv(f"{args.cs_path}/{args.gene}.indep.sushie.cs.tsv", sep="\t")
df_mega = pd.read_csv(f"{args.cs_path}/{args.gene}.normal.mega.cs.tsv", sep="\t")
df_meta = pd.read_csv(f"{args.cs_path}/{args.gene}.normal.meta.cs.tsv", sep="\t")

df_normal = df_normal[~df_normal.snp.isna()]
df_indep = df_indep[~df_indep.snp.isna()]
df_mega = df_mega[~df_mega.snp.isna()]
df_meta = df_meta[~df_meta.snp.isna()]

total_n = df_normal.shape[0] + df_indep.shape[0] + df_mega.shape[0] + df_meta.shape[0]

if total_n == 0:
    print("No credible set output for all the methods.")
    exit()
else:
    pd.DataFrame({"CS": ["Yes"]}).to_csv(args.cs_out, index=False)

