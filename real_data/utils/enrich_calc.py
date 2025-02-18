#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import statsmodels.api as sm
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
parser.add_argument("--weights-file", type=str)
parser.add_argument("--anno-file", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--n-l", type=int, default=10)
parser.add_argument("--exclude", type=str)
parser.add_argument("--enrich-out", type=str)

args = parser.parse_args()

df_normal = pd.read_csv(f"{args.weights_path}/sushie/alphas/{args.gene}.normal.sushie.alphas.tsv", sep="\t")
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

df_pips = pd.read_csv(args.weights_file, sep="\t")
df_pips.Value = correct_errors(df_pips.Value)
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
