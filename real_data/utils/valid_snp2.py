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

# def jaccard_index(set1, set2):
#     set1 = set(set1)
#     set2 = set(set2)
#     intersection = set1.intersection(set2)
#     union = set1.union(set2)
#     if len(union) == 0:
#         return 0  # Avoid division by zero when both sets are empty
#     return len(intersection) / len(union)

def get_cs_snps(df, col_name):
    df_sorted = df.sort_values(by=col_name, ascending=False)
    cumulative_sum = 0
    selected_snps = []

    # Iterate over the sorted DataFrame to include values until the cumulative sum reaches 0.95
    for index, row in df_sorted.iterrows():
        cumulative_sum += row[col_name]
        selected_snps.append(row["snp"])
        if cumulative_sum >= 0.95:
            break
    selected_snps_series = pd.Series(selected_snps)

    return selected_snps_series

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

df_w1 = pd.read_csv(args.file_w1, sep="\t")
df_w2 = pd.read_csv(args.file_w2, sep="\t")
df_w1.columns = df_w1.columns.str.lower()
df_w2.columns = df_w2.columns.str.lower()

all_res = []
for idx in range(1, 11):
    if f"alpha_l{idx}" in df_w1.columns:
        for jdx in range(1, 11):
            if f"alpha_l{jdx}" in df_w2.columns:
                main_snps = get_cs_snps(df_w1, f"alpha_l{idx}")
                valid_snps = get_cs_snps(df_w2, f"alpha_l{jdx}")
                cs_keep = main_snps.isin(valid_snps).any()
                df_s1 = df_w1[["snp", f"alpha_l{idx}"]]
                df_s2 = df_w2[["snp", f"alpha_l{jdx}"]]

                df_inner = df_s1.merge(df_s2, how="inner", on="snp")

                df_comp = df_s1.merge(df_s2, how="outer", on="snp").fillna(0)
                df_comp.columns = ["snp", "main_alpha", "valid_alpha"]

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


