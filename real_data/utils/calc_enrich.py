#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import statsmodels.api as sm
import warnings

def make_pip(alpha):
    pip = -np.expm1(np.sum(np.log1p(-alpha), axis=0))
    return pip

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--alpha-path", type=str)
parser.add_argument("--anno-file", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--n-l", type=int, default=10)
parser.add_argument("--exclude", type=str)
parser.add_argument("--enrich-out", type=str)

args = parser.parse_args()

df_normal = pd.read_csv(f"{args.alpha_path}/{args.gene}.normal.sushie.alphas.tsv", sep="\t")
df_anno = pd.read_csv(args.anno_file, sep="\t")
df_exc = pd.read_csv(args.exclude, sep="\t", header=None)
df_anno = df_anno[~df_anno.ANNOT.isin(df_exc[0].values)]

df_check = df_normal[["snp"]].copy()
ct = 0
for anno in df_anno.ANNOT.unique():
    tmp_check = df_anno[df_anno.ANNOT == anno].copy()
    df_check["y"] = 1 * (df_check["snp"].isin(tmp_check.SNP))
    if df_check.y.sum() == 0 or df_check.y.sum() == df_check.shape[0]:
        ct += 1

if ct == len(df_anno.ANNOT.unique()):
    print("All the SNPs have the same annotations 0 or 1.")
    exit()

df_indep = pd.read_csv(f"{args.alpha_path}/{args.gene}.indep.sushie.alphas.tsv", sep="\t")
df_mega = pd.read_csv(f"{args.alpha_path}/{args.gene}.normal.mega.alphas.tsv", sep="\t")
df_meta = pd.read_csv(f"{args.alpha_path}/{args.gene}.normal.meta.alphas.tsv", sep="\t")

# prepare pips, alphas, and annotation

res = []
# normal
normal_incs = np.sum(df_normal[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)

tmp_normal_pips = df_normal[["snp", "pip_all"]].copy()
tmp_normal_pips.columns = ["snp", "Value"]
tmp_normal_pips["Type"] = "pip_cov"
tmp_normal_pips = tmp_normal_pips[["snp", "Type", "Value"]]
res.append(tmp_normal_pips)
for ldx in range(args.n_l):
    if normal_incs[ldx] > 0:
        tmp_alpha_pips = df_normal[["snp", f"alpha_l{ldx + 1}"]].copy()
        tmp_alpha_pips.columns = ["snp", "Value"]
        tmp_alpha_pips["Type"] = f"alpha_l{ldx + 1}"
        tmp_alpha_pips = tmp_alpha_pips[["snp", "Type", "Value"]]
        res.append(tmp_alpha_pips)

# indep
tmp_indep_pips = df_indep[["snp", "pip_all"]].copy()
tmp_indep_pips.columns = ["snp", "Value"]
tmp_indep_pips["Type"] = "pip_indep"
tmp_indep_pips = tmp_indep_pips[["snp", "Type", "Value"]]
res.append(tmp_indep_pips)

# mega
tmp_mega_pips = df_mega[["snp", "pip_all"]].copy()
tmp_mega_pips.columns = ["snp", "Value"]
tmp_mega_pips["Type"] = "pip_mega"
tmp_mega_pips = tmp_mega_pips[["snp", "Type", "Value"]]
res.append(tmp_mega_pips)

# meta
meta_pips = df_meta[df_meta.ancestry == "ancestry_1"][["snp"]].copy()
for pop in df_meta.ancestry.unique():
    tmp_meta = df_meta[df_meta.ancestry == pop]
    tmp_meta_incs = np.sum(tmp_meta[[f"kept_l{idx}" for idx in range(1, args.n_l + 1)]] == 1)
    tmp_meta_pips = tmp_meta[["snp", "pip_all"]].copy()
    tmp_meta_pips.columns = ["snp", f"pip_{pop}"]
    meta_pips = meta_pips.merge(tmp_meta_pips, how="left", on="snp")
    tmp_meta_pips.columns = ["snp", "Value"]
    tmp_meta_pips["Type"] = f"pip_{pop}"
    tmp_meta_pips = tmp_meta_pips[["snp", "Type", "Value"]]
    res.append(tmp_meta_pips)

meta_pips["Type"] = "pip_meta"
n_pop = len(df_meta.ancestry.unique())
meta_pips["Value"] = make_pip(np.array(meta_pips.iloc[:, 1:(n_pop+1)].transpose()))
res.append(meta_pips[["snp", "Type", "Value"]])

df_pips = pd.concat(res)

df_y = df_normal[["snp"]].copy()
df_res = pd.DataFrame()
for anno in df_anno.ANNOT.unique():
    tmp_anno = df_anno[df_anno.ANNOT == anno].copy()
    df_y["y"] = 1 * (df_y["snp"].isin(tmp_anno.SNP))
    if df_y.y.sum() == 0 or df_y.y.sum() == df_y.shape[0]:
        continue
    else:
        for method_type in df_pips.Type.unique():
            tmp_pips = df_pips[df_pips.Type == method_type].merge(df_y, how="left", on="snp")
            X = sm.add_constant(np.array(tmp_pips.y))
            logit_model = sm.Logit(np.array(tmp_pips.Value), X)
            try:
                logit_model.fit(disp=False, warn_convergence=False)
            except Exception as e:
                print(f"Fit logistic regression with {method_type} on {anno} with regularization.")
                result = logit_model.fit_regularized(disp=False, warn_convergence=False)
                reg = 1
            else:
                result = logit_model.fit(disp=False, warn_convergence=False)
                reg = 0

            tmp_res = pd.DataFrame({"trait": [args.gene],
                                    "anno": anno,
                                    "method": [method_type],
                                    "est": [result.params[1].round(4)],
                                    "se": [result.bse[1].round(4)],
                                    "converged": [1*result.mle_retvals["converged"]],
                                    "reg": [reg],
                                    "n_snps": [tmp_pips.shape[0]],
                                    "n_anno": [df_y.y.sum()]})
            df_res = pd.concat([df_res, tmp_res])

df_res.to_csv(args.enrich_out, sep="\t", index=False)
