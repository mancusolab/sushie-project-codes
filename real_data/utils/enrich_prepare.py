#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import os
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp

def correct_errors(values, lb=0, ub=1):
    values = np.array(values)
    values = np.clip(values, lb, ub)
    return values

def make_pip(alpha):
    pip = -jnp.expm1(jnp.sum(jnp.log1p(-alpha), axis=0))
    return pip

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--weights-path", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--ans", type=int)
parser.add_argument("--n-l", type=int, default=10)
parser.add_argument("--out", type=str)

args = parser.parse_args()

df_normal = pd.read_csv(f"{args.weights_path}/sushie/alphas/{args.gene}.normal.sushie.alphas.tsv", sep="\t")
df_indep = pd.read_csv(f"{args.weights_path}/sushie/alphas/{args.gene}.indep.sushie.alphas.tsv", sep="\t")
df_mega = pd.read_csv(f"{args.weights_path}/sushie/alphas/{args.gene}.normal.mega.alphas.tsv", sep="\t")
df_meta = pd.read_csv(f"{args.weights_path}/sushie/alphas/{args.gene}.normal.meta.alphas.tsv", sep="\t")

df_normal_cs = pd.read_csv(f"{args.weights_path}/sushie/cs/{args.gene}.normal.sushie.cs.tsv", sep="\t")
df_indep_cs = pd.read_csv(f"{args.weights_path}/sushie/cs/{args.gene}.indep.sushie.cs.tsv", sep="\t")
df_mega_cs = pd.read_csv(f"{args.weights_path}/sushie/cs/{args.gene}.normal.mega.cs.tsv", sep="\t")
df_meta_cs = pd.read_csv(f"{args.weights_path}/sushie/cs/{args.gene}.normal.meta.cs.tsv", sep="\t")

# prepare pips, alphas, and annotation

res = []
# normal
if not df_normal_cs.snp.isna()[0]:
    tmp_normal_pips = df_normal[["snp", "pip_all"]].copy()
    tmp_normal_pips.columns = ["snp", "Value"]
    tmp_normal_pips["Type"] = "pip_cov"
    tmp_normal_pips = tmp_normal_pips[["snp", "Type", "Value"]]
    res.append(tmp_normal_pips)

    # # pip_cs
    # tmp_normal_pips = df_normal[["snp", "pip_cs"]].copy()
    # tmp_normal_pips.columns = ["snp", "Value"]
    # tmp_normal_pips["Type"] = "pip_cov_cs"
    # tmp_normal_pips = tmp_normal_pips[["snp", "Type", "Value"]]
    # res.append(tmp_normal_pips)

    normal_incs = np.sum(df_normal[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)
    for ldx in range(args.n_l):
        if normal_incs[ldx] > 0:
            tmp_alpha_pips = df_normal[["snp", f"alpha_l{ldx + 1}"]].copy()
            tmp_alpha_pips.columns = ["snp", "Value"]
            tmp_alpha_pips["Type"] = f"cov_l{ldx + 1}"
            res.append(tmp_alpha_pips[["snp", "Type", "Value"]])

# indep
if not df_indep_cs.snp.isna()[0]:
    tmp_indep_pips = df_indep[["snp", "pip_all"]].copy()
    tmp_indep_pips.columns = ["snp", "Value"]
    tmp_indep_pips["Type"] = "pip_indep"
    tmp_indep_pips = tmp_indep_pips[["snp", "Type", "Value"]]
    res.append(tmp_indep_pips)

    # tmp_indep_pips = df_indep[["snp", "pip_cs"]].copy()
    # tmp_indep_pips.columns = ["snp", "Value"]
    # tmp_indep_pips["Type"] = "pip_indep_cs"
    # tmp_indep_pips = tmp_indep_pips[["snp", "Type", "Value"]]
    # res.append(tmp_indep_pips)
    #
    # indep_incs = np.sum(df_indep[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)
    # for ldx in range(args.n_l):
    #     if indep_incs[ldx] > 0:
    #         tmp_alpha_pips = df_indep[["snp", f"alpha_l{ldx + 1}"]].copy()
    #         tmp_alpha_pips.columns = ["snp", "Value"]
    #         tmp_alpha_pips["Type"] = f"indep_l{ldx + 1}"
    #         res.append(tmp_alpha_pips[["snp", "Type", "Value"]])

# mega
if not df_mega_cs.snp.isna()[0]:
    tmp_mega_pips = df_mega[["snp", "pip_all"]].copy()
    tmp_mega_pips.columns = ["snp", "Value"]
    tmp_mega_pips["Type"] = "pip_mega"
    tmp_mega_pips = tmp_mega_pips[["snp", "Type", "Value"]]
    res.append(tmp_mega_pips)

    # tmp_mega_pips = df_mega[["snp", "pip_cs"]].copy()
    # tmp_mega_pips.columns = ["snp", "Value"]
    # tmp_mega_pips["Type"] = "pip_mega_cs"
    # tmp_mega_pips = tmp_mega_pips[["snp", "Type", "Value"]]
    # res.append(tmp_mega_pips)
    #
    # mega_incs = np.sum(df_mega[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)
    # for ldx in range(args.n_l):
    #     if mega_incs[ldx] > 0:
    #         tmp_alpha_pips = df_mega[["snp", f"alpha_l{ldx + 1}"]].copy()
    #         tmp_alpha_pips.columns = ["snp", "Value"]
    #         tmp_alpha_pips["Type"] = f"mega_l{ldx + 1}"
    #         res.append(tmp_alpha_pips[["snp", "Type", "Value"]])

# meta
if not df_meta_cs.snp.isna()[0]:
    meta_pips = df_meta[df_meta.ancestry == "ancestry_1"][["snp", "pip_all"]].copy()
    meta_pips.columns = ["snp", "Value"]
    meta_pips["Type"] = "pip_meta"
    res.append(meta_pips[["snp", "Type", "Value"]])

    # meta_pips = df_meta[df_meta.ancestry == "ancestry_1"][["snp", "pip_cs"]].copy()
    # meta_pips.columns = ["snp", "Value"]
    # meta_pips["Type"] = "pip_meta_cs"
    # res.append(meta_pips[["snp", "Type", "Value"]])
    #
    # meta_incs = np.sum(df_meta[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)
    # for ldx in range(args.n_l):
    #     if meta_incs[ldx] > 0:
    #         meta_pips = df_meta[df_meta.ancestry == "ancestry_1"][["snp"]].copy()
    #         meta_pips["Value"] = make_pip(df_meta[f"alpha_l{ldx + 1}"].values.reshape(df_normal.shape[0], args.ans).T)
    #         meta_pips["Type"] = f"meta_l{ldx + 1}"
    #         res.append(meta_pips[["snp", "Type", "Value"]])

susiex_paths = f"{args.weights_path}/susiex/weights/susiex.{args.gene}.weights.tsv"
mesusie_paths = f"{args.weights_path}/mesusie/weights/mesusie.{args.gene}.weights.tsv"
multisusie_paths = f"{args.weights_path}/multisusie/weights/multisusie.{args.gene}.weights.tsv"
# xmap_paths = f"{args.weights_path}/xmap/weights/xmap.{args.gene}.weights.tsv"

if os.path.exists(susiex_paths):
    tmp_df = pd.read_csv(susiex_paths, sep="\t")[["SNP", "pip_all"]].rename(columns={"SNP": "snp", "pip_all": "Value"}).copy()

    # pip_cs = jnp.empty((df_normal.shape[0], 0))
    # for idx in range(1, 11):
    #     if f"alpha_l{idx}" in tmp_df.columns:
    #         res_df = tmp_df[["SNP", f"alpha_l{idx}"]].copy()
    #         res_df.columns = ["snp", "Value"]
    #         res_df["Type"] = f"susiex_l{idx}"
    #         res_df = res_df[["snp", "Type", "Value"]]
    #         pip_cs = jnp.hstack((pip_cs, res_df.Value.values[:,jnp.newaxis]))
    #         res.append(res_df)
    # tmp_df2 = tmp_df[["SNP"]].rename(columns={"SNP": "snp"}).copy()
    # tmp_df2["Value"] = make_pip(pip_cs.T)
    tmp_df["Type"] = "susiex_cs"
    res.append(tmp_df[["snp", "Type", "Value"]])

if os.path.exists(mesusie_paths):
    tmp_df = pd.read_csv(mesusie_paths, sep="\t")[["snp", "pip_all"]].rename(columns={"pip_all": "Value"}).copy()
    # pip_cs = jnp.empty((df_normal.shape[0], 0))
    # for idx in range(1, 11):
    #     if f"alpha_l{idx}" in tmp_df.columns:
    #         res_df = tmp_df[["snp", f"alpha_l{idx}"]].copy()
    #         res_df.columns = ["snp", "Value"]
    #         res_df["Type"] = f"mesusie_l{idx}"
    #         res_df = res_df[["snp", "Type", "Value"]]
    #         pip_cs = jnp.hstack((pip_cs, res_df.Value.values[:,jnp.newaxis]))
    #         res.append(res_df)
    # tmp_df2 = tmp_df[["snp"]].copy()
    # tmp_df2["Value"] = make_pip(pip_cs.T)
    tmp_df["Type"] = "mesusie_cs"
    res.append(tmp_df[["snp", "Type", "Value"]])

# if os.path.exists(multisusie_paths):
#     tmp_df = pd.read_csv(multisusie_paths, sep="\t")
#     pip_cs = jnp.empty((df_normal.shape[0], 0))
#     for idx in range(1, 11):
#         if f"alpha_l{idx}" in tmp_df.columns:
#             res_df = tmp_df[["snp", f"alpha_l{idx}"]].copy()
#             res_df.columns = ["snp", "Value"]
#             res_df["Type"] = f"multisusie_l{idx}"
#             res_df = res_df[["snp", "Type", "Value"]]
#             pip_cs = jnp.hstack((pip_cs, res_df.Value.values[:,jnp.newaxis]))
#             res.append(res_df)
#     tmp_df2 = tmp_df[["snp"]].copy()
#     tmp_df2["Value"] = make_pip(pip_cs.T)
#     tmp_df2["Type"] = "multisusie_cs"
#     res.append(tmp_df2[["snp", "Type", "Value"]])


if len(res) == 0:
    quit()
df_pips = pd.concat(res)
df_pips.to_csv(args.out, sep="\t", index=False)
