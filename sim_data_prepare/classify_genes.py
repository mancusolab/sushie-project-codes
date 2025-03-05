#! /usr/bin/env python

import pandas as pd
import numpy as np
from pandas_plink import read_plink
import scipy.linalg as linalg
from scipy.stats import levene, ttest_ind, f_oneway

def _compute_ld(G):
    G = G.T
    # estimate LD for population from PLINK data
    n, p = [float(x) for x in G.shape]
    p_int = int(p)
    mafs = np.mean(G, axis=0) / 2
    G -= mafs * 2
    G /= np.std(G, axis=0)

    # regularize so that LD is PSD
    LD = np.dot(G.T, G) / n + np.eye(p_int) * 0.1

    # compute cholesky decomp for faster sampling/simulation
    # L = linalg.cholesky(LD, lower=True)
    L = linalg.cholesky(LD, lower=True, check_finite=False)
    mu = 2 * mafs

    adj_mu = linalg.solve_triangular(L, mu, lower=True, check_finite=False)

    return L, LD, adj_mu

def _allele_check(baseA0, baseA1, compareA0, compareA1):
    correct = np.array(
        ((baseA0 == compareA0) * 1) * ((baseA1 == compareA1) * 1), dtype=int
    )
    flipped = np.array(
        ((baseA0 == compareA1) * 1) * ((baseA1 == compareA0) * 1), dtype=int
    )
    correct_idx = np.where(correct == 1)[0]
    flipped_idx = np.where(flipped == 1)[0]
    wrong_idx = np.where((correct + flipped) == 0)[0]
    return correct_idx, flipped_idx, wrong_idx

df = pd.read_csv("/project/nmancuso_8/zeyunlu/projects/sushie-data-codes/sim_data_prepare/sim_gene_index_list.tsv",
                 sep="\t", header=None).sort_values(by=0)

path = "/project/nmancuso_8/data/sushie/analysis_results/sushie_sim_genes/new"

levene_stats = []
levene_p = []
anova_stats = []
anova_p = []
l2_diff = []
l2_diff_sq = []
n_snps = []
for idx in range(500):
    print(idx)
    bim = []
    fam = []
    bed = []
    n_pop = 2
    for jdx in range(n_pop):
        ans = "EUR" if jdx == 0 else "AFR"
        tmp_bim, tmp_fam, tmp_bed = read_plink(f"{path}/{ans}_1000G/{df[2].values[idx]}_{df[1].values[idx]}_geno", verbose=False)
        tmp_bim = tmp_bim[["chrom", "snp", "a0", "a1", "i"]].rename(
            columns={"i": f"bimIDX_{jdx}", "a0": f"a0_{jdx}", "a1": f"a1_{jdx}"})
        bim.append(tmp_bim)
        fam.append(tmp_fam)
        bed.append(tmp_bed)

    snps = pd.merge(bim[0], bim[1], how="inner", on=["chrom", "snp"])

    for jdx in range(n_pop - 2):
        snps = pd.merge(snps, bim[jdx + 2], how="inner", on=["chrom", "snp"])

    flip_jdx = []
    if n_pop > 1:
        for jdx in range(1, n_pop):
            # keep track of miss match alleles
            correct_jdx, tmp_flip_jdx, tmp_wrong_jdx = _allele_check(snps["a0_0"].values, snps["a1_0"].values,
                                                                     snps[f"a0_{jdx}"].values, snps[f"a1_{jdx}"].values)
            flip_jdx.append(tmp_flip_jdx)
            snps = snps.drop(index=tmp_wrong_jdx)
            snps = snps.drop(columns=[f"a0_{jdx}", f"a1_{jdx}"])
        snps = snps.rename(columns={"a0_0": "a0", "a1_0": "a1"})

    output_dic = {"LD": [], "L": [], "mu": []}

    all_ldsc = []
    LD_list = []
    for jdx in range(n_pop):
        bed[jdx] = bed[jdx][snps[f"bimIDX_{jdx}"].values, :].compute()

        # flip the mismatched allele
        if jdx > 0:
            bed[jdx][flip_jdx[jdx - 1]] = 2 - bed[jdx][flip_jdx[jdx - 1]]

        L, LD, mu = _compute_ld(bed[jdx])
        all_ldsc.append(np.sum(LD ** 2, axis=0))
        LD_list.append(LD)
    var_stats, var_p = levene(all_ldsc[0], all_ldsc[1])
    mean_stats, mean_p = f_oneway(all_ldsc[0], all_ldsc[1])
    l2_diff_val = np.abs(np.linalg.norm(LD_list[0]) - np.linalg.norm(LD_list[1]))
    l2_diff_sq_val = np.abs(np.linalg.norm(LD_list[0]**2) - np.linalg.norm(LD_list[1]**2))
    
    n_snps.append(bed[0].shape[0])
    levene_stats.append(var_stats)
    levene_p.append(var_p)
    anova_stats.append(mean_stats)
    anova_p.append(mean_p)
    l2_diff.append(l2_diff_val)
    l2_diff_sq.append(l2_diff_sq_val)

df["n_snps"] = np.array(n_snps)
df["levene_stats"] = np.array(levene_stats)
df["levene_p"] = np.array(levene_p)
df["anova_stats"] = np.array(anova_stats)
df["anova_p"] = np.array(anova_p)
df["l2_diff"] = np.array(l2_diff)
df["l2_diff_sq"] = np.array(l2_diff_sq)

df.to_csv("~/data/sushie/sim3/sim_gene_ld.tsv.gz", sep="\t", index=False)

