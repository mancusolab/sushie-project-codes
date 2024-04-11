#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse

def calc_bin_mean(df, strand, tss, tes, start, end, bin_edges, col):
    if strand == "+":
        df["dist"] = df.pos - tss
    else:
        df["dist"] = tes - df.pos

    df = df[df.dist <= end].copy()
    df = df[df.dist >= start].copy()
    df["bins"] = np.digitize(df.dist, bins=bin_edges)
    df = df.groupby("bins").agg(
        n_snps=("snp", "size"),
        mean=(col, "mean"),
    ).reset_index(names="bins")
    df["mean"] = df["mean"].round(8)
    df["Type"] = col
    return df


parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--alpha-path", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--strand", type=str)
parser.add_argument("--tss", type=int)
parser.add_argument("--tes", type=int)
parser.add_argument("--n-l", type=int, default=10)
parser.add_argument("--start", type=int, default=-5e5)
parser.add_argument("--end", type=int, default=5e5)
parser.add_argument("--step", type=int, default=500)
parser.add_argument("--tss-out", type=str)

args = parser.parse_args()

bin_edges = np.arange(args.start, (5e5+args.step), args.step)

# focus on sushie
df_normal = pd.read_csv(args.alpha_path, sep="\t")
df_incs = np.sum(df_normal[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)

# check pips
res = []
if np.sum(df_incs) > 0:
    tmp_pip = df_normal[["snp", "pos", "pip_all"]].copy()
    tmp_res = calc_bin_mean(tmp_pip, args.strand, args.tss, args.tes, args.start, args.end, bin_edges, "pip_all")
    res.append(tmp_res)
    for ldx in range(args.n_l):
        if df_incs[ldx] > 0:
            tmp_alpha = df_normal[["snp", "pos", f"alpha_l{ldx + 1}"]].copy()
            tmp_res = calc_bin_mean(tmp_alpha, args.strand, args.tss, args.tes, args.start, args.end, bin_edges, f"alpha_l{ldx + 1}")
            res.append(tmp_res)
    res = pd.concat(res)
    res["gene"] = args.gene
    res.to_csv(args.tss_out, sep="\t", index=False)
else:
    print("No credible sets are output.")
    exit()
