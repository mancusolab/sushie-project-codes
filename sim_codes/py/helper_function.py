import warnings
from pandas_plink import read_plink
import pandas as pd
import scipy.linalg as linalg
from sklearn import linear_model as lm
from numpy.linalg import multi_dot as mdot
from numpy_sugar.linalg import economic_qs
from scipy import stats
from jax.config import config
import numpy as np
from glimix_core.lmm import LMM
from jax import random

config.update("jax_enable_x64", True)

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp

def _compute_ld(G):
    G = G.T
    # estimate LD for population from PLINK data
    n, p = [float(x) for x in G.shape]
    p_int = int(p)
    mafs = jnp.mean(G, axis=0) / 2
    G -= mafs * 2
    G /= jnp.std(G, axis=0)

    # regularize so that LD is PSD
    original_LD = jnp.dot(G.T, G) / n
    LD = original_LD + jnp.eye(p_int) * 0.1

    # compute cholesky decomp for faster sampling/simulation
    # L = linalg.cholesky(LD, lower=True)
    L = linalg.cholesky(LD, lower=True, check_finite=False)
    mu = 2 * mafs
    adj_mu = linalg.solve_triangular(L, mu, lower=True, check_finite=False)
    # import pdb; pdb.set_trace()
    #
    # new_matrix = np.zeros_like(original_LD)
    # np.fill_diagonal(new_matrix, 1)
    # new_matrix[np.triu_indices_from(new_matrix, k=1)] = original_LD[np.triu_indices_from(original_LD, k=1)]
    # new_matrix[np.tril_indices_from(new_matrix, k=-1)] = new_matrix.T[np.tril_indices_from(new_matrix, k=-1)]
    return L, LD, adj_mu, original_LD

ambiguous_snps = ["AT", "TA", "CG", "GC"]

def _ambiguous_idx(baseA0, baseA1):
    if len(baseA0) != len(baseA1):
        raise ValueError("Base and compare alleles are not in the same dimension")

    ambiguous_idx = jnp.where((baseA0 + baseA1).isin(ambiguous_snps))[0]

    return ambiguous_idx

def _maf_idx(bim, bed, maf = 0.01):

    # calculate maf
    snp_maf = jnp.mean(bed, axis=1) / 2
    snp_maf = jnp.where(snp_maf > 0.5, 1 - snp_maf, snp_maf)

    sel_idx = jnp.where(snp_maf < maf)[0]

    if len(sel_idx) > 0:
        mask = jnp.ones(bed.shape[0], dtype=bool)
        mask = mask.at[sel_idx].set(False)
        bed = bed[mask, :]
        bim2 = bim.drop(sel_idx).reset_index(drop=True)
        bim2["i"] = bim2.index
        return bim2, bed
    else:
        return bim, bed

def _allele_check(baseA0, baseA1, compareA0, compareA1):

    if len(baseA0) != len(baseA1) or len(compareA0) != len(compareA1) or len(compareA0) != len(compareA1):
        raise ValueError("Base and compare alleles are not in the same dimension")

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


def _gen_ld(prefix_pop, remove_ambiguous=True):
    # read in plink data
    n_pop = len(prefix_pop)
    bim = []
    fam = []
    bed = []

    for idx in range(n_pop):
        tmp_bim, tmp_fam, tmp_bed = read_plink(prefix_pop[idx], verbose=False)
        tmp_bed = tmp_bed.compute()
        tmp_bim, tmp_bed = _maf_idx(tmp_bim, tmp_bed)

        tmp_bim = tmp_bim[["chrom", "snp", "pos", "a0", "a1", "i"]].rename(
            columns={"i": f"bimIDX_{idx}", "a0": f"a0_{idx}", "a1": f"a1_{idx}"})

        bim.append(tmp_bim)
        fam.append(tmp_fam)
        bed.append(tmp_bed)

    snps = pd.merge(bim[0], bim[1], how="inner", on=["chrom", "snp"])

    for idx in range(n_pop - 2):
        snps = pd.merge(snps, bim[idx + 2], how="inner", on=["chrom", "snp"])

    snps = snps.reset_index(drop=True)
    if remove_ambiguous:
        # we just need to drop ambiguous index for the first ancestry
        # and then it will automatically drop for others in later codes
        ambig_idx = _ambiguous_idx(snps["a0_0"].values, snps["a1_0"].values)
        snps = snps.drop(index=ambig_idx).reset_index(drop=True)

    # remove wrong alleles
    if n_pop > 1:
        for idx in range(1, n_pop):
            _, _, tmp_wrong_idx = _allele_check(snps["a0_0"].values, snps["a1_0"].values,
                                                snps[f"a0_{idx}"].values, snps[f"a1_{idx}"].values)
            snps = snps.drop(index=tmp_wrong_idx).reset_index(drop=True)

    flip_idx = []
    if n_pop > 1:
        for idx in range(1, n_pop):
            correct_idx, tmp_flip_idx, _ = _allele_check(snps["a0_0"].values, snps["a1_0"].values,
                                                         snps[f"a0_{idx}"].values, snps[f"a1_{idx}"].values)
            flip_idx.append(tmp_flip_idx)
            snps = snps.drop(columns=[f"a0_{idx}", f"a1_{idx}"])
        snps = snps.rename(columns={"a0_0": "a0", "a1_0": "a1"})

    output_dic = {"LD": [], "L": [], "mu": [], "trueLD": []}
    for idx in range(n_pop):
        bed[idx] = bed[idx][snps[f"bimIDX_{idx}"].values, :]
        # flip the mismatched allele
        if idx > 0 and len(flip_idx[idx - 1]) != 0:
            bed[idx][flip_idx[idx - 1]] = 2 - bed[idx][flip_idx[idx - 1], :]

        L, LD, mu, trueLD = _compute_ld(bed[idx])
        output_dic["L"].append(L)
        output_dic["LD"].append(LD)
        output_dic["mu"].append(mu)
        output_dic["trueLD"].append(trueLD)

    return output_dic, snps

def simulation_sushie(rng_key, output_dic, N, L, L3, h2g, rho):
    mu = output_dic["mu"]
    L_cho = output_dic["L"]

    p, _ = L_cho[0].shape
    n_pop = len(N)

    X = []
    b_var = []

    for idx in range(n_pop):
        rng_key, x_key = random.split(rng_key, 2)
        tmp_X = L_cho[idx].dot((random.normal(x_key, shape=(N[idx], p)) + mu[idx]).T).T
        tmp_X -= jnp.mean(tmp_X, axis=0)
        tmp_X /= jnp.std(tmp_X, axis=0)
        X.append(tmp_X)
        b_var.append(h2g[idx] / (L + L3))

    b_covar = jnp.diag(jnp.array(b_var))
    ct = 0
    for row in range(1, n_pop):
        for col in range(n_pop):
            if col < row:
                _cov = jnp.sqrt(b_var[row] * b_var[col])
                b_covar = b_covar.at[row, col].set(rho[ct] * _cov)
                b_covar = b_covar.at[col, row].set(rho[ct] * _cov)
                ct += 1

    ct = 0
    # Make sure we have signficant region
    while ct < 100:
        ct = ct + 1
        rng_key, b_key = random.split(rng_key, 2)

        if jnp.sum(b_covar) != 0:
            bvec = random.multivariate_normal(b_key, jnp.zeros((n_pop,)), b_covar, shape=(L,))
        else:
            bvec = jnp.zeros((L, n_pop))

        b_indep = jnp.zeros((n_pop, L3))
        for idx in range(n_pop):
            rng_key, b_indep_key = random.split(rng_key, 2)
            b_indep = b_indep.at[idx, :].set(jnp.sqrt(b_var[idx]) * random.normal(b_indep_key, shape=(L3,)))

        rng_key, gamma_key = random.split(rng_key, 2)
        gamma = random.choice(gamma_key, p, shape=(L + n_pop * L3,), replace=False)

        bvec_all = []
        for idx in range(n_pop):
            start_pos = L + L3 * idx
            end_pos = start_pos + L3 * (idx + 1)
            tmp_bvec = jnp.zeros(p).at[gamma[0:L]].set(jnp.transpose(bvec)[idx]).at[gamma[start_pos:end_pos]].set(
                b_indep[idx])
            bvec_all.append(tmp_bvec)

        g = []
        s2e = []
        y = []
        lrt = []
        for idx in range(n_pop):
            tmp_g = X[idx] @ bvec_all[idx]
            tmp_s2g = jnp.var(tmp_g)

            if h2g[idx] != 0:
                tmp_s2e = ((1 / h2g[idx]) - 1) * tmp_s2g
            else:
                tmp_s2e = 1

            rng_key, y_key = random.split(rng_key, 2)
            tmp_y = tmp_g + jnp.sqrt(tmp_s2e) * random.normal(y_key, shape=(N[idx],))

            _, _, _, tmp_lrt, _ = estimate_her(X[idx], tmp_y)

            tmp_y -= jnp.mean(tmp_y)
            tmp_y /= jnp.std(tmp_y)
            g.append(tmp_g)
            s2e.append(tmp_s2e)
            y.append(tmp_y)
            lrt.append(tmp_lrt)

        # 2.7055
        keep_greater = True
        # for idx in range(n_pop):
        #     keep_greater = False if lrt[idx] < 2.7055 else keep_greater

        if keep_greater:
            break

    return rng_key, X, y, bvec, bvec_all, s2e, gamma, b_covar

def make_pip(alpha):
    pip = -jnp.expm1(jnp.sum(jnp.log1p(-alpha), axis=0))
    return pip

def ols(X: jnp.ndarray, y: jnp.ndarray):
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


def sim_trait(g, h2g):
    n = len(g)

    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ((1.0 / h2g) - 1)
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e
    else:
        e = np.random.normal(0, 1, n)
        y = e

    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()

    return y, y_std


def regress(Z, pheno):
    betas = []
    ses = []
    pvals = []
    for snp in Z.T:
        beta, inter, rval, pval, se = stats.linregress(snp, pheno)
        betas.append(beta)
        ses.append(se)
        pvals.append(pval)

    res = pd.DataFrame({"beta": betas, "se": ses, "pval": pvals})

    return res


def sim_gwas(L, ngwas, b_qtls, var_explained, x_key):
    Z_gwas = sim_geno(L, int(ngwas), x_key)

    # var_explained should only reflect that due to genetics
    gwas_expr = np.dot(Z_gwas, b_qtls)
    if var_explained > 0:
        alpha = np.random.normal(loc=0, scale=1)
    else:
        alpha = 0

    y, y_std = sim_trait(gwas_expr * alpha, var_explained)

    gwas = regress(Z_gwas, y)

    # correct alpha for original SD of y
    alpha /= y_std

    return (gwas, alpha)


def sim_geno(L, n, x_key):
    p, p = L.shape

    Z = L.dot(random.normal(x_key, shape=(n, p)).T).T
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0)

    return Z


def compute_twas(gwas, coef, LD):
    Z = gwas.beta.values / gwas.se.values

    # score and variance
    score = np.einsum("ji,i->j", coef.T, Z)
    within_var = np.einsum("jk,ki,ij->j", coef.T, LD, coef)

    z_twas = np.full(len(score), np.nan)
    # p_twas = np.ones(len(score))
    for idx in range(len(score)):
        if within_var[idx] > 0:
            z_twas[idx] = score[idx] / np.sqrt(within_var[idx])

    return z_twas
