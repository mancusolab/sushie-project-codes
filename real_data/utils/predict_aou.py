#! /usr/bin/env python

import argparse
from typing import Tuple

import jax.numpy as jnp
import numpy as np
import pandas as pd
from glimix_core.lmm import LMM
from jax import Array
import jax.scipy as jsp
from jax.typing import ArrayLike
from numpy.linalg import multi_dot as mdot
from numpy_sugar.linalg import economic_qs
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm
from jax import random
from sushie.infer import infer_sushie
import os

def fit_lasso(rng_key, Z, y, h2g):
    return _fit_sparse_penalized_model(rng_key, Z, y, h2g, lm.Lasso)


def fit_enet(rng_key, Z, y, h2g):
    return _fit_sparse_penalized_model(rng_key, Z, y, h2g, lm.ElasticNet)


def fit_ridge(rng_key, Z, y, h2g):
    n, p = Z.shape
    lambda_r = (1 - h2g) / (h2g / p)

    model = lm.Ridge(alpha=np.array(lambda_r))
    model.fit(np.array(Z), np.array(y))
    coef, r2, logl = _get_model_info(model, Z, y)

    return coef, r2, logl


def _fit_sparse_penalized_model(rng_key, Z, y, h2g, model_cls=lm.Lasso):
    np.random.seed(rng_key)

    if model_cls not in [lm.Lasso, lm.ElasticNet]:
        raise ValueError("penalized model must be either Lasso or ElasticNet")

    n, p = Z.shape

    def _gen_e():
        e = np.random.normal(size=n)
        return np.linalg.norm(Z.T.dot(e), np.inf)

    # PLINK-style LASSO
    lambda_max = np.linalg.norm(Z.T.dot(y), np.inf) / float(n)

    min_tmp = np.median([_gen_e() for _ in range(1000)])
    sige = np.sqrt(1.0 - h2g + (1.0 / float(n)))
    lambda_min = (sige / n) * min_tmp

    # 100 values spaced logarithmically from lambda-min to lambda-max
    alphas = np.exp(np.linspace(np.log(lambda_min), np.log(lambda_max), 100))

    # fit solution using coordinate descent, updating with consecutively smaller penalties
    model = model_cls(fit_intercept=True, warm_start=True)
    for penalty in reversed(alphas):
        model.set_params(alpha=penalty)
        model.fit(Z, y)

    coef, r2, logl = _get_model_info(model, Z, y)

    return coef, r2, logl


def _get_model_info(model, Z, y):
    n, p = Z.shape
    coef = model.coef_

    r2 = model.score(Z, y)
    ystar = model.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl

def _allele_check(baseA0, baseA1, compareA0, compareA1):
    correct = jnp.array(
        ((baseA0 == compareA0) * 1) * ((baseA1 == compareA1) * 1), dtype=int
    )
    flipped = jnp.array(
        ((baseA0 == compareA1) * 1) * ((baseA1 == compareA0) * 1), dtype=int
    )
    correct_idx = jnp.where(correct == 1)[0]
    flipped_idx = jnp.where(flipped == 1)[0]
    wrong_idx = jnp.where((correct + flipped) == 0)[0]

    return correct_idx, flipped_idx, wrong_idx

def ols(X: ArrayLike, y: ArrayLike) -> Tuple[Array, Array, Array]:
    X_inter = jnp.append(jnp.ones((X.shape[0], 1)), X, axis=1)
    y = jnp.reshape(y, (len(y), -1))
    q_matrix, r_matrix = jnp.linalg.qr(X_inter, mode="reduced")
    qty = q_matrix.T @ y
    beta = jsp.linalg.solve_triangular(r_matrix, qty)
    df = q_matrix.shape[0] - q_matrix.shape[1]
    residual = y - q_matrix @ qty
    rss = jnp.sum(residual ** 2, axis=0)
    sigma = jnp.sqrt(jnp.sum(residual ** 2, axis=0) / df)
    se = (
            jnp.sqrt(
                jnp.diag(
                    jsp.linalg.cho_solve((r_matrix, False), jnp.eye(r_matrix.shape[0]))
                )
            )[:, jnp.newaxis]
            @ sigma[jnp.newaxis, :]
    )
    t_scores = beta / se
    p_value = jnp.array(2 * stats.t.sf(abs(t_scores), df=df))

    r_sq = 1 - rss / jnp.sum((y - jnp.mean(y, axis=0)) ** 2, axis=0)
    adj_r = 1 - (1 - r_sq) * (q_matrix.shape[0] - 1) / df

    return residual, adj_r, p_value

def estimate_her(X, y):
    n, p = X.shape

    covar = jnp.ones(n)

    GRM = jnp.dot(X, X.T) / p
    GRM = GRM / jnp.diag(GRM).mean()
    QS = economic_qs(GRM)
    method = LMM(y, covar, QS, restricted=True)
    method.fit(verbose=False)

    g = method.scale * (1 - method.delta)
    e = method.scale * method.delta
    v = jnp.var(method.mean())
    h2g_w_v = g / (v + g + e)
    h2g_wo_v = g / (g + e)
    alt_lk = method.lml()
    method.delta = 1
    method.fix("delta")
    method.fit(verbose=False)
    null_lk = method.lml()
    lrt_stats = -2 * (null_lk - alt_lk)
    p_value = stats.chi2.sf(lrt_stats, 1) / 2

    return g, h2g_w_v, h2g_wo_v, lrt_stats, p_value

def regress_covar(
    X: ArrayLike, y: ArrayLike, covar: ArrayLike, no_regress: bool
) -> Tuple[Array, Array]:

    y, _, _ = ols(covar, y)
    if not no_regress:
        X, _, _ = ols(covar, X)

    return X, y


parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--pheno", type=str, nargs="+")
parser.add_argument("--covar", type=str, nargs="+")
parser.add_argument("--plink", type=str, nargs="+")
parser.add_argument("--trait", type=str)
parser.add_argument("--w_files", type=str)
parser.add_argument("--seed", type=int)
parser.add_argument("--study", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

n_pop = len(args.plink)
sushie_cs_path = f"{args.w_files}/sushie/cs/{args.trait}.normal.sushie.cs.tsv"
indep_cs_path = f"{args.w_files}/sushie/cs/{args.trait}.indep.sushie.cs.tsv"
mega_cs_path = f"{args.w_files}/sushie/cs/{args.trait}.normal.mega.cs.tsv"
meta_cs_path = f"{args.w_files}/sushie/cs/{args.trait}.normal.meta.cs.tsv"
mesusie_cs_path = f"{args.w_files}/mesusie/cs/mesusie.{args.trait}.cs.tsv"

na_sushie = pd.isna(pd.read_csv(sushie_cs_path, sep="\t")["snp"].iloc[0])
na_indep = pd.isna(pd.read_csv(indep_cs_path, sep="\t")["snp"].iloc[0])
na_meta = pd.isna(pd.read_csv(meta_cs_path, sep="\t")["snp"].iloc[0])
na_mega = pd.isna(pd.read_csv(mega_cs_path, sep="\t")["snp"].iloc[0])
if os.path.exists(mesusie_cs_path):
    na_mesusie = False
else:
    na_mesusie = True

if na_sushie and na_indep and na_mega and na_meta and na_mesusie:
    print("No CS for all methods")
    exit()

bim = []
fam = []
bed = []

# prepare geno
for idx in range(n_pop):
    tmp_bim, tmp_fam, tmp_bed = read_plink(args.plink[idx], verbose=False)
    tmp_bim = tmp_bim[["chrom", "snp", "pos", "a0", "a1", "i"]].rename(
        columns={"i": f"bimIDX_{idx}", "a0": f"a0_{idx}", "a1": f"a1_{idx}"})
    bim.append(tmp_bim)
    fam.append(tmp_fam)
    bed.append(tmp_bed.compute())

snps = pd.merge(bim[0], bim[1], how="inner", on=["chrom", "snp", "pos"])
chr_num = snps["chrom"][0]

for idx in range(n_pop - 2):
    snps = pd.merge(snps, bim[idx + 2], how="inner", on=["chrom", "snp", "pos"])

flip_idx = []
for idx in range(1, n_pop):
    correct_idx, tmp_flip_idx, _ = _allele_check(snps["a0_0"].values, snps["a1_0"].values,
                                                 snps[f"a0_{idx}"].values, snps[f"a1_{idx}"].values)
    flip_idx.append(tmp_flip_idx)
    snps = snps.drop(columns=[f"a0_{idx}", f"a1_{idx}"])
snps = snps.rename(columns={"a0_0": "a0", "a1_0": "a1"})

for idx in range(n_pop):
    bed[idx] = bed[idx][snps[f"bimIDX_{idx}"].values, :]
    # flip the mismatched allele
    if idx > 0 and len(flip_idx[idx - 1]) != 0:
        bed[idx][flip_idx[idx - 1]] = 2 - bed[idx][flip_idx[idx - 1], :]

pheno = []
geno = []
for idx in range(n_pop):
    tmp_pheno = pd.read_csv(args.pheno[idx], sep="\t", header=None)
    tmp_covar = pd.read_csv(args.covar[idx], sep="\t", header=None)
    tmp_pheno = fam[idx][["iid"]].merge(tmp_pheno, how="inner", left_on="iid", right_on=0).drop(columns=0)
    tmp_pheno.columns = ["iid", "pheno"]
    tmp_covar = fam[idx][["iid"]].merge(tmp_covar, how="inner", left_on="iid", right_on=0).drop(columns=0)
    bed[idx] = bed[idx].T - jnp.mean(bed[idx].T, axis=0)
    bed[idx] /= jnp.std(bed[idx], axis=0)
    new_bed, new_pheno = regress_covar(bed[idx], tmp_pheno.pheno.values, jnp.array(tmp_covar.drop(columns="iid")), False)
    new_bed -= jnp.mean(new_bed, axis=0)
    new_pheno -= jnp.mean(new_pheno)
    new_bed /= jnp.std(new_bed, axis=0)
    new_pheno /= jnp.std(new_pheno)
    pheno.append(jnp.squeeze(new_pheno))
    geno.append(new_bed)

comb_bed = jnp.concatenate(geno)
comb_pheno = jnp.concatenate(pheno)

comb_bed -= jnp.mean(comb_bed, axis=0)
comb_pheno -= jnp.mean(comb_pheno)
comb_bed /= jnp.std(comb_bed, axis=0)
comb_pheno /= jnp.std(comb_pheno)

_, h2g, _, _, _ = estimate_her(comb_bed, comb_pheno)

enet1, _, _ = fit_enet(args.seed, comb_bed, comb_pheno, h2g)

lasso1, _, _ = fit_lasso(args.seed, comb_bed, comb_pheno, h2g)

ridge1, _, _ = fit_ridge(args.seed, comb_bed, comb_pheno, h2g)

snps["enet_weight"] = enet1
snps["lasso_weight"] = lasso1
snps["gblup_weight"] = ridge1
snps["rsid"] = snps.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
snps = snps[["rsid", "a1", "enet_weight", "lasso_weight", "gblup_weight"]]

df_tmp = pd.read_csv(f"{args.w_files}/sushie/weights/{args.trait}.normal.sushie.weights.tsv", sep="\t")
if n_pop == 2:
    df_tmp["ancestry3_sushie_weight"] = 0

if na_sushie:
    df_tmp["ancestry1_sushie_weight"] = 0
    df_tmp["ancestry2_sushie_weight"] = 0
    df_tmp["ancestry3_sushie_weight"] = 0

sel_col = ["chrom", "pos", "a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight", "ancestry3_sushie_weight", "snp"]
sel_col2 = ["rsid", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight",
           "ancestry3_sushie_weight", "snp"]
rename_col = ["rsid", "a1", "sushie_pop1", "sushie_pop2", "sushie_pop3", "snp"]
df_tmp = df_tmp[sel_col]
df_tmp["rsid"] = df_tmp.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp = df_tmp[sel_col2]
df_tmp.columns = rename_col

df_tmp_mega = pd.read_csv(f"{args.w_files}/sushie/weights/{args.trait}.normal.mega.weights.tsv", sep="\t")
if na_mega:
    df_tmp_mega["mega_weight"] = 0

df_tmp_mega = df_tmp_mega[["chrom", "pos","a0", "a1", "mega_weight"]]
df_tmp_mega["rsid"] = df_tmp_mega.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp_mega = df_tmp_mega[["rsid", "a1", "mega_weight"]]

df_tmp_meta = pd.read_csv(f"{args.w_files}/sushie/weights/{args.trait}.normal.meta.weights.tsv", sep="\t")
if n_pop == 2:
    df_tmp_meta["ancestry3_single_weight"] = 0
if na_meta:
    df_tmp_meta["ancestry1_single_weight"] = 0
    df_tmp_meta["ancestry2_single_weight"] = 0
    df_tmp_meta["ancestry3_single_weight"] = 0
df_tmp_meta = df_tmp_meta[["chrom", "pos", "a0", "a1", "ancestry1_single_weight", "ancestry2_single_weight", "ancestry3_single_weight"]]
df_tmp_meta["rsid"] = df_tmp_meta.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp_meta = df_tmp_meta[["rsid", "a1", "ancestry1_single_weight", "ancestry2_single_weight", "ancestry3_single_weight"]]
df_tmp_meta.columns = ["rsid", "a1", "meta_pop1", "meta_pop2", "meta_pop3"]

df_tmp_indep = pd.read_csv(f"{args.w_files}/sushie/weights/{args.trait}.indep.sushie.weights.tsv", sep="\t")
if n_pop == 2:
    df_tmp_indep["ancestry3_sushie_weight"] = 0

if na_indep:
    df_tmp_indep["ancestry1_sushie_weight"] = 0
    df_tmp_indep["ancestry2_sushie_weight"] = 0
    df_tmp_indep["ancestry3_sushie_weight"] = 0

df_tmp_indep = df_tmp_indep[["chrom", "pos", "a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight", "ancestry3_sushie_weight"]]
df_tmp_indep["rsid"] = df_tmp_indep.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp_indep = df_tmp_indep[["rsid", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight", "ancestry3_sushie_weight"]]
df_tmp_indep.columns = ["rsid", "a1", "indep_pop1", "indep_pop2", "indep_pop3"]

mesusie_paths = f"{args.w_files}/mesusie/weights/mesusie.{args.trait}.weights.tsv"

if os.path.exists(mesusie_paths):
    df_tmp_mesusie = pd.read_csv(mesusie_paths, sep="\t")
    if n_pop == 2:
        df_tmp_mesusie["pop3_weights"] = 0
    df_tmp_mesusie = df_tmp_mesusie[["snp", "pop1_weights", "pop2_weights", "pop3_weights"]]
    df_tmp_mesusie.columns = ["snp", "mesusie_pop1", "mesusie_pop2", "mesusie_pop3"]
else:
    df_tmp_mesusie = df_tmp[["snp"]].copy()
    df_tmp_mesusie["mesusie_pop1"] = 0
    df_tmp_mesusie["mesusie_pop2"] = 0
    df_tmp_mesusie["mesusie_pop3"] = 0

# multisusie_paths = f"{args.w_files}/multisusie/weights/multisusie.{args.trait}.weights.tsv"
# if os.path.exists(multisusie_paths):
#     df_tmp_multisusie = pd.read_csv(multisusie_paths, sep="\t")
#     if n_pop == 2:
#         df_tmp_multisusie["pop3_weights"] = 0
#     df_tmp_multisusie = df_tmp_multisusie[["snp", "pop1_weights", "pop2_weights", "pop3_weights"]]
#     df_tmp_multisusie.columns = ["snp", "multisusie_pop1", "multisusie_pop2", "multisusie_pop3"]
# else:
#     df_tmp_multisusie = df_tmp[["snp"]].copy()
#     df_tmp_multisusie["multisusie_pop1"] = 0
#     df_tmp_multisusie["multisusie_pop2"] = 0
#     df_tmp_multisusie["multisusie_pop3"] = 0

df_res = df_tmp.merge(df_tmp_indep, how="inner", on=["rsid", "a1"]) \
    .merge(df_tmp_meta, how="inner", on=["rsid", "a1"]) \
        .merge(df_tmp_mega, how="inner", on=["rsid", "a1"]) \
        .merge(df_tmp_mesusie, how="inner", on="snp").merge(snps, how="inner", on=["rsid", "a1"])
df_res = df_res.drop(columns=["snp"])

df_res.to_csv(f"{args.out}/weights/chr{chr_num}/{args.study}.chr{chr_num}.{args.trait}.weights.tsv", sep="\t", index=False)
df_res[["rsid"]].to_csv(f"{args.out}/snps/chr{chr_num}/{args.study}.chr{chr_num}.{args.trait}.snps.tsv", sep="\t", index=False, header=None)
