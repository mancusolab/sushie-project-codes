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
parser.add_argument("--group", type=int)
parser.add_argument("--study", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

n_pop = len(args.plink)
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
df_tmp = pd.read_csv(f"{args.w_files}/{args.trait}.group{args.group}.normal.sushie.weights.tsv", sep="\t")
snps = snps[snps.snp.isin(df_tmp.snp)]
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
snps = snps[["snp", "a1", "enet_weight", "lasso_weight", "gblup_weight"]]

df_tmp = pd.read_csv(f"{args.w_files}/{args.trait}.group{args.group}.normal.sushie.weights.tsv", sep="\t")
if n_pop == 2:
    df_tmp["ancestry3_sushie_weight"] = 0
df_tmp = df_tmp[["snp", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight", "ancestry3_sushie_weight"]]
df_tmp.columns = ["snp", "a1", "sushie_pop1", "sushie_pop2", "sushie_pop3"]

df_tmp_mega = pd.read_csv(f"{args.w_files}/{args.trait}.group{args.group}.normal.mega.weights.tsv", sep="\t")
df_tmp_mega = df_tmp_mega[["snp", "a1", "mega_weight"]]

df_res = df_tmp.merge(df_tmp_mega, how="inner", on=["snp", "a1"]).merge(snps, how="inner", on=["snp", "a1"])

# $TMPDIR/${ID}.EUR.test.group${idx}.geno
ans_list = ["EUR", "AFR", "HIS"] if n_pop == 3 else ["EUR", "AFR"]

for idx in range(n_pop):
    tmp_bim, tmp_fam, tmp_bed = read_plink(f"{args.w_files}/{args.trait}.{ans_list[idx]}.test.group{args.group}.geno", verbose=False)
    tmp_bim = tmp_bim[["chrom", "snp", "a0", "a1", "i"]]
    tmp_bim = tmp_bim.merge(df_res[["snp"]], how="right", on="snp")
    tmp_bed = tmp_bed.compute()[tmp_bim.i.values, :].T
    tmp_res = df_res[["snp", "a1", f"sushie_pop{idx+1}", "mega_weight", "enet_weight", "lasso_weight", "gblup_weight"]]
    df_wk = tmp_bim[["snp", "a0", "a1"]].merge(tmp_res, how="inner", on=["snp"])
    comp = ((df_wk.a1_x == df_wk.a1_y) * 1).replace(0, -1)
    all_weights = jnp.einsum("ij,i->ij", jnp.array(df_wk.iloc[:,4:9]), jnp.array(comp.astype(int)))
    imp_pheno = jnp.einsum("ij,jk->ik", tmp_bed, all_weights)
    tmp_res2 = pd.concat([tmp_fam[["iid"]], pd.DataFrame(imp_pheno)], axis=1)
    tmp_res2.columns = ["IID", "sushie", "mega", "enet", "lasso", "ridge"]
    tmp_res2.to_csv(f"{args.out}.{ans_list[idx]}.value.sscore", sep="\t", index=False)


