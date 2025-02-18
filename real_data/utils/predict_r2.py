#! /usr/bin/env python

import argparse
import jax.numpy as jnp
import pandas as pd
import os

def calc_r2(X: jnp.ndarray, y: jnp.ndarray):
    n_samples, n_features = X.shape
    X_inter = jnp.append(jnp.ones((n_samples, 1)), X, axis=1)
    n_features += 1
    y = jnp.reshape(y, (len(y), -1))
    XtX_inv = jnp.linalg.inv(X_inter.T @ X_inter)
    betas = XtX_inv @ X_inter.T @ y
    residual = y - X_inter @ betas
    rss = jnp.sum(residual ** 2, axis=0)
    r_sq = 1 - rss / jnp.sum((y - jnp.mean(y, axis=0)) ** 2, axis=0)
    adj_r = 1 - (1 - r_sq) * (n_samples - 1) / (n_samples - n_features)

    return r_sq, adj_r


parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--ans", type=int)
parser.add_argument("--trait", type=str)
parser.add_argument("--wk_folder", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

n_pop = args.ans

ans_list = ["EUR", "AFR", "HIS"] if n_pop == 3 else ["EUR", "AFR"]

import statsmodels.api as sm
import statsmodels.tools.tools as sm_tools
import numpy as np

all_df_r2 = []
for idx in range(n_pop):
    pop_value = []
    pop_cross = []
    for jdx in range(1, 6):
        df_tmp = pd.read_csv(f"{args.wk_folder}/{args.trait}.{ans_list[idx]}.test.group{jdx}.pheno", sep="\t", header=None)
        df_tmp["ans"] = ans_list[idx]
        df_tmp.columns = ["IID", "pheno", "ans"]
        path1=f"{args.wk_folder}/imp.pheno.{args.trait}.group{jdx}.{ans_list[idx]}.value.sscore"
        path2=f"{args.wk_folder}/imp.pheno.{args.trait}.group{jdx}.{ans_list[idx]}.cross.sscore"
        if os.path.exists(path1) and os.path.exists(path2):
            df_tmp2 = pd.read_csv(path1, sep="\t")
            df_tmp3 = pd.read_csv(path2, sep="\t")

            df_tmp_value = df_tmp[["IID", "ans", "pheno"]].merge(df_tmp2, how="inner", on=["IID"])
            df_tmp_cross = df_tmp[["IID", "ans", "pheno"]].merge(df_tmp3, how="inner", on=["IID"])
            pop_value.append(df_tmp_value)
            pop_cross.append(df_tmp_cross)
    
    if len(pop_value) == 0 or len(pop_cross) == 0:
        quit()
    
    df_pop_value = pd.concat(pop_value, axis=0)
    df_pop_cross = pd.concat(pop_cross, axis=0)

    X1 = sm_tools.add_constant(np.array(df_pop_value.iloc[:, 2]))
    X2 = sm_tools.add_constant(np.array(df_pop_cross.iloc[:, 2]))
    r2_1 = []
    r2_2 = []
    r2p_1 = []
    r2p_2 = []
    for col in range(3, 11):
        model1 = sm.OLS(np.array(df_pop_value.iloc[:, col]), X1).fit()
        model2 = sm.OLS(np.array(df_pop_cross.iloc[:, col]), X2).fit()
        r2_1.append(model1.rsquared)
        r2_2.append(model2.rsquared)
        r2p_1.append(model1.f_pvalue)
        r2p_2.append(model2.f_pvalue)

    df_r2 = pd.DataFrame(np.array(r2_1)[jnp.newaxis,:])
    df_r2.columns = ["sushie", "indep", "meta", "susie", "mesusie", "enet", "lasso", "ridge"]
    df_r2["trait"] = args.trait
    df_r2["group"] = "match"
    df_r2["type"] = "r2"

    df_r2p = pd.DataFrame(np.array(r2p_1)[jnp.newaxis, :])
    df_r2p.columns = ["sushie", "indep", "meta", "susie", "mesusie", "enet", "lasso", "ridge"]
    df_r2p["trait"] = args.trait
    df_r2p["group"] = "match"
    df_r2p["type"] = "r2_p"

    df_r2_2 = pd.DataFrame(np.array(r2_2)[jnp.newaxis, :])
    df_r2_2.columns = ["sushie", "indep", "meta", "susie", "mesusie", "enet", "lasso", "ridge"]
    df_r2_2["trait"] = args.trait
    df_r2_2["group"] = "cross"
    df_r2_2["type"] = "r2"

    df_r2p_2 = pd.DataFrame(np.array(r2p_2)[jnp.newaxis, :])
    df_r2p_2.columns = ["sushie", "indep", "meta", "susie", "mesusie", "enet", "lasso", "ridge"]
    df_r2p_2["trait"] = args.trait
    df_r2p_2["group"] = "cross"
    df_r2p_2["type"] = "r2_p"

    df_r2 = pd.concat([df_r2, df_r2p, df_r2_2, df_r2p_2], axis=0)
    df_r2["pop"] = ans_list[idx]
    all_df_r2.append(df_r2)

pd.concat(all_df_r2, axis=0).to_csv(args.out, sep="\t", index=False)

