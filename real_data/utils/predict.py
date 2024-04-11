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

def rint(y_val: ArrayLike) -> Array:
    n_pt = y_val.shape[0]
    r_y = stats.rankdata(y_val)
    q_y = stats.norm.ppf(r_y / (n_pt + 1))

    return q_y

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
parser.add_argument("--out", type=str)
parser.add_argument("--out_r2", type=str)

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
    bed.append(tmp_bed)

snps = pd.merge(bim[0], bim[1], how="inner", on=["chrom", "snp", "pos"])

for idx in range(n_pop - 2):
    snps = pd.merge(snps, bim[idx + 2], how="inner", on=["chrom", "snp", "pos"])

flip_idx = []
if n_pop > 1:
    for idx in range(1, n_pop):
        # keep track of miss match alleles
        correct_idx, tmp_flip_idx, tmp_wrong_idx = _allele_check(snps["a0_0"].values, snps["a1_0"].values,
                                                                 snps[f"a0_{idx}"].values, snps[f"a1_{idx}"].values)
        flip_idx.append(tmp_flip_idx)
        snps = snps.drop(index=tmp_wrong_idx)
        snps = snps.drop(columns=[f"a0_{idx}", f"a1_{idx}"])
    snps = snps.rename(columns={"a0_0": "a0", "a1_0": "a1"})

for idx in range(n_pop):
    bed[idx] = bed[idx][snps[f"bimIDX_{idx}"].values, :].compute()
    # flip the mismatched allele
    if idx > 0:
        bed[idx][flip_idx[idx - 1]] = 2 - bed[idx][flip_idx[idx - 1]]

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

_, h2g, _, _, _ = estimate_her(comb_bed, comb_pheno)

enet1, _, _ = fit_enet(args.seed, comb_bed, comb_pheno, h2g)

lasso1, _, _ = fit_lasso(args.seed, comb_bed, comb_pheno, h2g)

ridge1, _, _ = fit_ridge(args.seed, comb_bed, comb_pheno, h2g)

snps["enet_weight"] = enet1
snps["lasso_weight"] = lasso1
snps["gblup_weight"] = ridge1
snps["rsid"] = snps.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
snps = snps[["rsid", "a1", "enet_weight", "lasso_weight", "gblup_weight"]]

df_tmp = pd.read_csv(f"{args.w_files}.normal.sushie.weights.tsv", sep="\t")
if n_pop == 3:
    sel_col = ["chrom", "pos","a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight", "ancestry3_sushie_weight"]
    sel_col2 = ["rsid", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight",
               "ancestry3_sushie_weight"]
    rename_col = ["rsid", "a1", "EUR_weight", "AFR_weight", "HIS_weight"]
else:
    sel_col = ["chrom", "pos", "a0", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight"]
    sel_col2 = ["rsid", "a1", "ancestry1_sushie_weight", "ancestry2_sushie_weight"]
    rename_col = ["rsid", "a1", "EUR_weight", "AFR_weight"]

df_tmp = df_tmp[sel_col]
df_tmp["rsid"] = df_tmp.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp = df_tmp[sel_col2]
df_tmp.columns = rename_col

df_tmp2 = pd.read_csv(f"{args.w_files}.normal.mega.weights.tsv", sep="\t")
df_tmp2 = df_tmp2[["chrom", "pos","a0", "a1", "mega_weight"]]
df_tmp2["rsid"] = df_tmp2.apply(lambda x: f"chr{x['chrom']}:{x['pos']}:{x['a1']}:{x['a0']}", axis=1)
df_tmp2 = df_tmp2[["rsid", "a1", "mega_weight"]]

df_res = df_tmp.merge(df_tmp2, how="inner", on=["rsid","a1"]).merge(snps, how="inner", on=["rsid","a1"])
df_res.to_csv(args.out, sep="\t", index=False)

rng_key = random.PRNGKey(args.seed)
cv_num = 5
geno_split = []
pheno_split = []
# shuffle the data first
for idx in range(n_pop):
    tmp_n = geno[idx].shape[0]
    rng_key, c_key = random.split(rng_key, 2)
    shuffled_index = random.choice(c_key, tmp_n, (tmp_n,), replace=False)
    geno[idx] = geno[idx][shuffled_index]
    pheno[idx] = pheno[idx][shuffled_index]
    geno_split.append(jnp.array_split(geno[idx], cv_num))
    pheno_split.append(jnp.array_split(pheno[idx], cv_num))

ori_pheno = []
est_pheno = []
est_pheno_cross = []
for cv in range(cv_num):
    train_geno = []
    train_pheno = []
    valid_geno = []
    valid_pheno = []
    train_index = jnp.delete(jnp.arange(cv_num), cv).tolist()

    # make the training and test for each population separately
    # because sample size may be different
    for idx in range(n_pop):
        valid_geno_tmp = geno_split[idx][cv]
        valid_geno_tmp -= jnp.mean(valid_geno_tmp, axis=0)
        valid_geno_tmp /= jnp.std(valid_geno_tmp, axis=0)
        valid_geno.append(valid_geno_tmp)
        train_geno_tmp = jnp.concatenate([geno_split[idx][jdx] for jdx in train_index])
        train_geno_tmp -= jnp.mean(train_geno_tmp, axis=0)
        train_geno_tmp /= jnp.std(train_geno_tmp, axis=0)
        train_geno.append(train_geno_tmp)

        valid_pheno_tmp = pheno_split[idx][cv]
        valid_pheno_tmp -= jnp.mean(valid_pheno_tmp)
        valid_pheno_tmp /= jnp.std(valid_pheno_tmp)
        valid_pheno.append(valid_pheno_tmp)

        train_pheno_tmp = jnp.concatenate([pheno_split[idx][jdx] for jdx in train_index])
        train_pheno_tmp -= jnp.mean(train_pheno_tmp)
        train_pheno_tmp /= jnp.std(train_pheno_tmp)
        train_pheno.append(train_pheno_tmp)

    ori_pheno.append(valid_pheno)
    sushie = infer_sushie(train_geno, train_pheno)
    prior_rho = [0] if n_pop == 2 else [0, 0, 0]
    indep = infer_sushie(train_geno, train_pheno, rho=prior_rho, no_update=True)

    # susie
    rb_X = jnp.concatenate(train_geno)
    rb_y = jnp.concatenate(train_pheno)

    susie = infer_sushie([rb_X], [rb_y])

    _, h2g, _, _, _ = estimate_her(rb_X, rb_y)
    enet1, _, _ = fit_enet(args.seed, rb_X, rb_y, h2g)
    lasso1, _, _ = fit_lasso(args.seed, rb_X, rb_y, h2g)
    ridge1, _, _ = fit_ridge(args.seed, rb_X, rb_y, h2g)

    sushie_weight = jnp.sum(sushie.posteriors.post_mean, axis=0)
    indep_weight = jnp.sum(indep.posteriors.post_mean, axis=0)
    susie_weight = jnp.sum(susie.posteriors.post_mean, axis=0)

    # meta
    meta_weights = []
    for idx in range(n_pop):
        tmp_meta = infer_sushie([train_geno[idx]], [train_pheno[idx]])
        tmp_weight = jnp.sum(tmp_meta.posteriors.post_mean, axis=0)
        meta_weights.append(tmp_weight)

    est_y = []
    est_y_cross = []
    for idx in range(n_pop):
        tmp_all_weights = jnp.concatenate((sushie_weight[:, idx][:, jnp.newaxis],
                                           indep_weight[:, idx][:, jnp.newaxis],
                                           meta_weights[idx], susie_weight,
                                           enet1[:, jnp.newaxis], lasso1[:, jnp.newaxis],
                                           ridge1[:, jnp.newaxis]), axis=1)
        est_y1 = jnp.einsum("ij,jk->ik", valid_geno[idx], tmp_all_weights)
        est_y.append(est_y1)
        if n_pop == 2:
            if idx == 0:
                jdx = 1
            else:
                jdx = 0
        else:
            if idx == 0:
                jdx = 1
            elif idx == 1:
                jdx = 2
            else:
                jdx = 0
        est_y2 = jnp.einsum("ij,jk->ik", valid_geno[idx], sushie_weight[:, jdx][:, jnp.newaxis])
        est_y_cross.append(est_y2)

    est_pheno.append(est_y)
    est_pheno_cross.append(est_y_cross)

all_ori = []
all_est = []
all_est_cross = []
for idx in range(n_pop):
    tmp_ori = jnp.concatenate([inner_list[idx] for inner_list in ori_pheno])
    all_ori.append(tmp_ori)
    tmp_est = jnp.concatenate([inner_list[idx] for inner_list in est_pheno])
    all_est.append(tmp_est)
    tmp_est_cross = jnp.concatenate([inner_list[idx] for inner_list in est_pheno_cross])
    all_est_cross.append(tmp_est_cross)

df_ori = jnp.concatenate(all_ori)
df_est = jnp.concatenate(all_est)
df_est_cross = jnp.concatenate(all_est_cross)
r2, _ = calc_r2(df_ori[:, np.newaxis], df_est)
r2_cross, _ = calc_r2(df_ori[:, np.newaxis], df_est_cross)

df_r2 = pd.DataFrame(r2).T
df_r2_cross = pd.DataFrame(r2_cross)
df_res = pd.concat([df_r2, df_r2_cross], axis=1)
df_res.columns = ["sushie", "indep", "meta", "susie", "enet", "lasso", "ridge", "cross"]
df_res["type"] = "r2"
df_res["trait"] = args.trait
df_res.to_csv(args.out_r2, sep="\t", index=False)
