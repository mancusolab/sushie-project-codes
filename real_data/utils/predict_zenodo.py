#! /usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--trait", type=str)
parser.add_argument("--w_files", type=str)
parser.add_argument("--study", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

sushie_cs_path = f"{args.w_files}/sushie/cs/{args.trait}.normal.sushie.cs.tsv"
na_sushie = pd.isna(pd.read_csv(sushie_cs_path, sep="\t")["snp"].iloc[0])

if na_sushie:
    print("No CS")
    exit()

df_tmp = pd.read_csv(f"{args.w_files}/sushie/weights/{args.trait}.normal.sushie.weights.tsv", sep="\t")
if args.study == "genoa.mrna":
    sel_col = ["chrom", "snp", "pos", "a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight"]
    rename_col = ["chrom", "snp", "pos", "a0", "a1", "EUR_weight", "AFR_weight"]
else:
    sel_col = ["chrom", "snp", "pos", "a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight",
               "ancestry3_sushie_weight"]
    rename_col = ["chrom", "snp", "pos", "a0", "a1", "EUR_weight", "AFR_weight", "HIS_weight"]

df_tmp = df_tmp[sel_col]
df_tmp.columns = rename_col
df_tmp["trait"] = args.trait
df_tmp["study"] = args.study

df_tmp.to_csv(args.out, sep="\t", index=False)
