#! /usr/bin/env python
import warnings
import argparse
from pandas_plink import read_plink
import pandas as pd
import scipy.linalg as linalg
from jax.config import config
import numpy as np
import MultiSuSiE

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

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--ss_file", type=str)
parser.add_argument("--geno_file", type=str)
parser.add_argument("--wgt_file", type=str)
parser.add_argument("--rm_amb", type=str)
parser.add_argument("--trait", type=str)
parser.add_argument("--tmp_out", type=str)
parser.add_argument("--out", type=str)

args = parser.parse_args()

ss_path = args.ss_file.split(":")
plink_path = args.geno_file.split(":")
df_wgt = pd.read_csv(f"{args.wgt_file}", sep="\t")

rm_ambigous = True if args.rm_amb == "True" else False

output_dic, snps, sample_size = _gen_ld(plink_path, df_wgt.snp.values, remove_ambiguous=rm_ambigous)

n_pop = len(plink_path)

df_ss = pd.read_csv(f"{ss_path[0]}", sep="\t")

tmp_ss_list = []
for idx in range(n_pop):
    df_ss = pd.read_csv(f"{ss_path[idx]}", sep="\t")
    df_ss = df_ss[["ID", "A1", "BETA", "SE", "T_STAT", "P"]]
    df_z = snps.merge(df_ss, how="left", left_on="snp", right_on="ID").reset_index(drop=True)
    comp = ((df_z.A1 == df_z.a1) * 1).replace(0, -1)
    df_z["T_STAT"] = df_z["T_STAT"] * comp
    df_z["BETA"] = df_z["BETA"] * comp
    df_z = df_z[["chrom", "snp", "pos_x", "a0", "a1", "T_STAT", "BETA", "SE", "P"]].rename(columns={"pos_x": "pos"})
    new_dfz = df_z[~np.isnan(df_z.T_STAT)].copy()
    df_z = df_z[~np.isnan(df_z.T_STAT)].reset_index(drop=True)
    tmp_ss_list.append(df_z)
    output_dic["trueLD"][idx] = output_dic["trueLD"][idx][new_dfz.index.values, :][:, new_dfz.index.values]

common_snps = tmp_ss_list[0][["snp"]].merge(tmp_ss_list[1][["snp"]], how="inner", on="snp")
common_snps = common_snps.merge(tmp_ss_list[2][["snp"]], how="inner", on="snp")

ss_list = []
new_ld_list = []
for idx in range(n_pop):
    new_z1 = tmp_ss_list[idx][tmp_ss_list[idx].snp.isin(common_snps.snp.values)].copy()
    new_z2 = tmp_ss_list[idx][tmp_ss_list[idx].snp.isin(common_snps.snp.values)].reset_index(drop=True)
    ss_list.append(new_z2)
    new_z2.to_csv(f"{args.tmp_out}.inss.ans{idx}.tsv", index=False, sep="\t")
    tmp_ld = output_dic["trueLD"][idx][new_z1.index.values, :][:, new_z1.index.values]
    new_ld_list.append(tmp_ld)
    ld_df = pd.DataFrame(tmp_ld)
    ld_df.columns = new_z2.snp.values
    ld_df.to_csv(f"{args.tmp_out}.inld.ans{idx}.tsv", index=False, sep="\t", float_format="%.12f")

# compute ldscore for xmap and mesusie
ldsc1 = np.sum(new_ld_list[0] ** 2, axis=0)
ldsc2 = np.sum(new_ld_list[1] ** 2, axis=0)
ldsc3 = np.sum(new_ld_list[2] ** 2, axis=0)
ldsc12 = np.sum(new_ld_list[0] * new_ld_list[1], axis=0)
ldsc13 = np.sum(new_ld_list[0] * new_ld_list[2], axis=0)
ldsc23 = np.sum(new_ld_list[1] * new_ld_list[2], axis=0)
pd.DataFrame(jnp.concatenate([ldsc1[:, jnp.newaxis], ldsc2[:, jnp.newaxis], ldsc3[:, jnp.newaxis],
                             ldsc12[:, jnp.newaxis], ldsc13[:, jnp.newaxis], ldsc23[:, jnp.newaxis]], axis=1)). \
    to_csv(f"{args.tmp_out}.inldsc.tsv", index=False, sep="\t", header=None)

df_met = pd.DataFrame({"N": sample_size,
                    "BP": [jnp.min(df_z.pos.values), jnp.max(df_z.pos.values), 0]})

df_met.to_csv(f"{args.tmp_out}.md.tsv", index=False, sep="\t", header=None)

# multisusie1_converged = True
# try:
#     multisusie1 = MultiSuSiE.multisusie_rss(z_list=[np.array(ss_list[0].T_STAT.values), np.array(ss_list[1].T_STAT.values), np.array(ss_list[2].T_STAT.values)],
#                                             R_list=new_ld_list[0:3], rho=np.array([[1, 0.1, 0.1], [0.1, 1, 0.1], [0.1, 0.1, 1]]),
#                                             population_sizes=sample_size, L=10, max_iter=500, tol=0.0001, min_abs_corr=0.5)
# except Exception as e:
#     multisusie1_converged = False
#
# # cs
# multisusie1_cs = []
# for idx in range(10):
#     if multisusie1_converged and multisusie1.sets[3][idx]:
#         df_tmp = pd.DataFrame({"CSIndex": idx + 1, "SNPIndex": multisusie1.sets[0][idx]})
#         df_tmp["alpha"] = multisusie1.alpha.T[multisusie1.sets[0][idx],idx]
#         multisusie1_cs.append(df_tmp)
#
# if len(multisusie1_cs) != 0:
#     multisusie1_cs = pd.concat(multisusie1_cs, axis=0)
#     df_tmp1 = ss_list[0][["snp"]].merge(snps[["chrom", "snp", "pos_x", "a0", "a1"]], how="left", on="snp").reset_index(
#         names=["SNPIndex"])
#     df_cs = df_tmp1.merge(multisusie1_cs, how="right", on="SNPIndex")
#     df_cs = df_cs.rename(columns={"pos_x": "pos"})
#     df_cs["trait"] = args.trait
#     df_cs = df_cs[["chrom", "snp", "pos", "a0", "a1", "CSIndex", "SNPIndex", "alpha", "trait"]]
#     df_cs["method"] = "multisusie"
#     df_cs.to_csv(f"{args.out}/cs/multisusie.{args.trait}.cs.tsv", index=False, sep="\t")
#
#     multisusie1_coef_pop1 = multisusie1.coef[0][:, np.newaxis]
#     multisusie1_coef_pop2 = multisusie1.coef[1][:, np.newaxis]
#     multisusie1_coef_pop3 = multisusie1.coef[2][:, np.newaxis]
#     df_weights = snps[snps.snp.isin(common_snps.snp.values)][["chrom", "snp", "pos_x", "a0", "a1"]].copy().rename(columns={"pos_x": "pos"})
#     df_weights["pop1_weights"] = multisusie1_coef_pop1
#     df_weights["pop2_weights"] = multisusie1_coef_pop2
#     df_weights["pop3_weights"] = multisusie1_coef_pop3
#     df_weights["pip"] = multisusie1.pip
#     for idx in range(10):
#         if pd.Series([idx + 1]).isin(multisusie1_cs.CSIndex.values).any():
#             df_weights[f"alpha_l{idx + 1}"] = multisusie1.alpha.T[:, idx]
#
#     df_weights["trait"] = args.trait
#     df_weights["method"] = "multisusie"
#     df_weights.to_csv(f"{args.out}/weights/multisusie.{args.trait}.weights.tsv", index=False, sep="\t", float_format="%.12f")



