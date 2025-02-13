#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from numpy.linalg import norm
import copy

def make_pip(alpha):
    pip = -np.expm1(np.sum(np.log1p(-alpha), axis=0))
    return pip

def calc_p(z_stat):
    cdf = 0.5 * (1 + np.math.erf(z_stat / np.sqrt(2)))
    p_value = (1 - cdf)
    return p_value

def cosine_similarity(arr1, arr2):
    dot_product = np.dot(arr1, arr2)
    norm_arr1 = norm(arr1)
    norm_arr2 = norm(arr2)
    if norm_arr1 == 0 or norm_arr2 == 0:
        return 0
    cosine_sim = dot_product / (norm_arr1 * norm_arr2)
    return cosine_sim

# def weighted_jaccard_numpy(A, B):
#     min_sum = np.sum(np.minimum(A, B))
#     max_sum = np.sum(np.maximum(A, B))
#     return min_sum / max_sum if max_sum != 0 else 0
#
# def jaccard_index(set1, set2):
#     set1 = set(set1)
#     set2 = set(set2)
#     intersection = set1.intersection(set2)
#     union = set1.union(set2)
#     if len(union) == 0:
#         return 0  # Avoid division by zero when both sets are empty
#     return len(intersection) / len(union)


def meta_ancestry(df_wk, n_ans, n_range):
    df_snps = df_wk[df_wk.ancestry == "ancestry_1"][["snp"]].copy()
    for kdx in range(1, n_range + 1):
        df_snps1 = df_wk[df_wk.ancestry == "ancestry_1"][["snp"]].copy()
        df_snps2 = df_wk[df_wk.ancestry == "ancestry_1"][["snp"]].copy()
        df_snps3 = df_wk[df_wk.ancestry == "ancestry_1"][["snp"]].copy()
        for jdx in range(1, n_ans + 1):
            df_tmp = df_wk[df_wk.ancestry == f"ancestry_{jdx}"].copy()
            tmp_keep = df_tmp[f"kept_l{kdx}"].values[0] == 1
            if tmp_keep:
                df_snps1 = df_snps1.merge(df_tmp[["snp", f"alpha_l{kdx}"]], how="inner", on="snp")
                df_snps2 = df_snps2.merge(df_tmp[["snp", f"in_cs_l{kdx}"]], how="inner", on="snp")
                df_snps3 = df_snps3.merge(df_tmp[["snp", f"kept_l{kdx}"]], how="inner", on="snp")

        df_snps1[f"alpha_l{kdx}"] = make_pip(np.array(df_snps1.iloc[:, 1:(n_ans + 1)].transpose()))
        df_snps2[f"in_cs_l{kdx}"] = (np.sum(df_snps2.iloc[:, 1:(n_ans + 1)], axis=1) > 0) * 1
        df_snps3[f"kept_l{kdx}"] = (np.sum(df_snps3.iloc[:, 1:(n_ans + 1)], axis=1) > 0) * 1

        df_snps1 = df_snps1[["snp", f"alpha_l{kdx}"]]\
            .merge(df_snps2[["snp", f"in_cs_l{kdx}"]], how="inner", on="snp")\
            .merge(df_snps3[["snp", f"kept_l{kdx}"]], how="inner", on="snp")

        df_snps = df_snps.merge(df_snps1, how="inner", on="snp")

    return df_snps

parser = argparse.ArgumentParser(description="extract phenotype")

parser.add_argument("--file_w1", type=str)
parser.add_argument("--file_w2", type=str)
parser.add_argument("--gene1", type=str)
parser.add_argument("--gene2", type=str)
parser.add_argument("--seed", type=int)
parser.add_argument("--rep", default=100, type=int)
parser.add_argument("--ans", default=0, type=int)
parser.add_argument("--out", type=str)

args = parser.parse_args()

np.random.seed(args.seed)

df_o1 = pd.read_csv(args.file_w1, sep="\t")
df_o2 = pd.read_csv(args.file_w2, sep="\t")

# to see if we work on meta
if args.ans != 0:
    df_w1 = meta_ancestry(df_o1, args.ans, 10)
    df_w2 = meta_ancestry(df_o2, args.ans, 10)
else:
    df_w1 = df_o1
    df_w2 = df_o2

all_res = []
for idx in range(1, 11):
    w1_keep = df_w1[f"kept_l{idx}"][0] == 1
    if w1_keep:
        for jdx in range(1, 11):
            w2_keep = df_w2[f"kept_l{jdx}"][0] == 1
            if w2_keep:
                w2_tmp = df_w2[df_w2[f"in_cs_l{jdx}"] == 1]
                cs_keep = df_w1[df_w1[f"in_cs_l{idx}"] == 1].snp.isin(w2_tmp.snp).any()

                df_s1 = df_w1[["snp", f"alpha_l{idx}", f"in_cs_l{idx}"]]
                df_s2 = df_w2[["snp", f"alpha_l{jdx}", f"in_cs_l{jdx}"]]

                main_snps = df_s1[df_s1[f"in_cs_l{idx}"] == 1].snp.unique()
                valid_snps = df_s2[df_s2[f"in_cs_l{jdx}"] == 1].snp.unique()

                df_inner = df_s1.merge(df_s2, how="inner", on="snp")

                df_comp = df_s1.merge(df_s2, how="outer", on="snp").fillna(0)
                df_comp.columns = ["snp", "main_alpha", "main_cs", "valid_alpha", "valid_cs"]

                # w_jac = weighted_jaccard_numpy(df_comp.main_alpha.values, df_comp.valid_alpha.values)
                # n_jac = jaccard_index(main_snps, valid_snps)
                cos = cosine_similarity(df_comp.main_alpha.values, df_comp.valid_alpha.values)

                # null_w_jac = []
                # null_n_jac = []
                null_cos = []
                for ldx in range(args.rep):
                    test_pip = copy.deepcopy(df_comp.valid_alpha.values)
                    np.random.shuffle(test_pip)
                    # tmp_w_jac = weighted_jaccard_numpy(df_comp.main_alpha.values, test_pip)
                    tmp_cos = cosine_similarity(df_comp.main_alpha.values, test_pip)
                    # tmp_n_jac = jaccard_index(main_snps, np.random.choice(df_w2.snp.values, size=len(valid_snps), replace=False))
                    # null_w_jac.append(tmp_w_jac)
                    null_cos.append(tmp_cos)
                    # null_n_jac.append(tmp_n_jac)

                # m_w_jac = np.mean(null_w_jac)
                # sd_w_jac = np.std(null_w_jac)
                # z_w_jac = (w_jac - m_w_jac) / sd_w_jac
                # p_w_jac = calc_p(z_w_jac)
                #
                # m_n_jac = np.mean(null_n_jac)
                # sd_n_jac = np.std(null_n_jac)
                # z_n_jac = (n_jac - m_n_jac) / sd_n_jac
                # p_n_jac = calc_p(z_n_jac)

                m_cos = np.mean(null_cos)
                sd_cos = np.std(null_cos)
                z_cos = (cos - m_cos) / sd_cos
                p_cos = calc_p(z_cos)

                res = pd.DataFrame({"main_gene": [args.gene1],
                                    "valid_gene": [args.gene2],
                                    "main_cs_index": [idx],
                                    "valid_cs_index": [jdx],
                                    "match": [cs_keep * 1],
                                    "n_main": [df_w1.shape[0]],
                                    "n_valid": [df_w2.shape[0]],
                                    "n_overlap": [df_inner.shape[0]],
                                    "n_union": [df_comp.shape[0]],
                                    "n_cs_main": [len(main_snps)],
                                    "n_cs_valid": [len(valid_snps)],
                                    # "w_jac": [w_jac],
                                    # "mean_w_jac": [m_w_jac],
                                    # "se_w_jac": [sd_w_jac],
                                    # "z_w_jac": [z_w_jac],
                                    # "p_w_jac": [p_w_jac],
                                    # "n_jac": [n_jac],
                                    # "mean_n_jac": [m_n_jac],
                                    # "se_n_jac": [sd_n_jac],
                                    # "z_n_jac": [z_n_jac],
                                    # "p_n_jac": [p_n_jac],
                                    "cos": [cos],
                                    "mean_cos": [m_cos],
                                    "se_cos": [sd_cos],
                                    "z_cos": [z_cos],
                                    "p_cos": [p_cos],
                                    })
                all_res.append(res)

if len(all_res) != 0:
    pd.concat(all_res).to_csv(args.out, sep="\t", index=False)
else:
    print("no CS for 10 effects")


