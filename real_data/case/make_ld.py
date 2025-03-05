#! /usr/bin/env python
import warnings
from pandas_plink import read_plink
import pandas as pd
import scipy.linalg as linalg
from jax.config import config

config.update("jax_enable_x64", True)

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

    L = linalg.cholesky(LD, lower=True, check_finite=False)
    mu = 2 * mafs
    adj_mu = linalg.solve_triangular(L, mu, lower=True, check_finite=False)

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

def _gen_ld(prefix_pop, cand_snps, remove_ambiguous=True):
    # read in plink data
    n_pop = len(prefix_pop)
    bim = []
    fam = []
    bed = []
    sample_size = []
    for idx in range(n_pop):
        tmp_bim, tmp_fam, tmp_bed = read_plink(prefix_pop[idx], verbose=False)
        tmp_bed = tmp_bed.compute()
        tmp_bim = tmp_bim[tmp_bim["snp"].isin(cand_snps)]
        tmp_bed = tmp_bed[tmp_bim["i"].values, :]
        tmp_bim = tmp_bim.reset_index(drop=True).reset_index()
        tmp_bim["i"] = tmp_bim["index"]
        tmp_bim = tmp_bim.drop(columns=["index"])
        tmp_bim, tmp_bed = _maf_idx(tmp_bim, tmp_bed)

        tmp_bim = tmp_bim[["chrom", "snp", "pos", "a0", "a1", "i"]].rename(
            columns={"i": f"bimIDX_{idx}", "a0": f"a0_{idx}", "a1": f"a1_{idx}"})
        sample_size.append(tmp_fam.shape[0])
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

    return output_dic, snps, sample_size

# annoying warnings if on Mac ARM M1
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import jax.numpy as jnp

# gene_name = "ENSG00000106608_URGCP"
gene_name = "ENSG00000203875_SNHG5"

plink_path = [f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.EUR.geno",
              f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.AFR.geno",
              f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.HIS.geno"]

df_wgt = pd.read_csv(f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.normal.sushie.weights.tsv", sep="\t")

rm_ambigous = True

# snp_name = "rs2528382"
snp_name = "rs1059307"
output_dic, snps, sample_size = _gen_ld(plink_path, df_wgt.snp.values, remove_ambiguous=rm_ambigous)
snp_index = snps[snps["snp"] == snp_name].index.values
snp_ld = snps[["snp"]]
for idx in range(3):
    snp_ld["ldr2"] = output_dic["trueLD"][idx][:, snp_index] ** 2
    snp_ld.to_csv(f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.ans{idx}.{snp_name}.ldr2.tsv", sep="\t", index=False)

import pdb; pdb.set_trace()
# Define the window size
window_size = 25000

# Sort the SNPs by position
snps = snps.sort_values(by="pos_x").reset_index(drop=True)

# Group the SNPs by intervals of 5000 bp
snps["group"] = (snps["pos_x"] - snps["pos_x"].min()) // window_size

# Get the number of unique groups
num_groups = snps["group"].nunique()

# Initialize a new LD matrix for the grouped SNPs

# Create the new LD matrix by averaging correlations for each group
for kdx in range(3):
    new_ld_matrix = []# Diagonal is 1
    for idx in range(num_groups):
        for jdx in range(idx, num_groups):
            if idx == jdx:
                new_ld_matrix.append(pd.DataFrame({"pos1": [jnp.min(snps[snps["group"] == idx].pos_x.values)], "pos2": [jnp.min(snps[snps["group"] == jdx].pos_x.values)], "counts": [1]}))
                continue
            # Get the indices of SNPs in groups i and j
            group_i_indices = snps[snps["group"] == idx].index
            group_j_indices = snps[snps["group"] == jdx].index
            # Extract the submatrix from the original LD matrix
            submatrix = output_dic["trueLD"][kdx][jnp.array(group_i_indices),:][:,jnp.array(group_j_indices)]

            # Calculate the average correlation for this group pair
            avg_correlation = jnp.mean(submatrix**2)
            new_ld_matrix.append(pd.DataFrame({"pos1": [jnp.min(snps[snps["group"] == idx].pos_x.values)], "pos2": [jnp.min(snps[snps["group"] == jdx].pos_x.values)], "counts": [avg_correlation]}))
    pd.concat(new_ld_matrix).to_csv(f"/Users/zeyunlu/Documents/github/data/sushie_results/case/{gene_name}/{gene_name}.{kdx}.trueLD.tsv", sep="\t", index=False)
