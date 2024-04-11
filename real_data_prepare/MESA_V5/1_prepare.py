#!/usr/bin/env python3

import pandas as pd
import numpy as np
import qtl.norm

path = "/project/nmancuso_8/data/MESA/dbgap/WGS/rna_seq_v2"
scretch_path = "/scratch1/zeyunlu/sushie"

pt1 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010512.v1.p1.c2.TOPMed_MESA_RNA_Seq_Expression_Sample_Attributes.HMB-NPU.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "expr_id",
                                           "AGE_AT_COLLECTION": "age",
                                           "HISTOLOGICAL_TYPE": "type",
                                           "COLLECTION_VISIT": "visit",
                                           "ASSAY_LAB": "assay_lab"})
pt1 = pt1[["expr_id", "type", "visit", "assay_lab", "age"]]

pt2 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010512.v1.p1.c1.TOPMed_MESA_RNA_Seq_Expression_Sample_Attributes.HMB.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "expr_id",
                                           "AGE_AT_COLLECTION": "age",
                                           "HISTOLOGICAL_TYPE": "type",
                                           "COLLECTION_VISIT": "visit",
                                           "ASSAY_LAB": "assay_lab"})
pt2 = pt2[["expr_id", "type", "visit", "assay_lab", "age"]]
df_pt = pd.concat([pt1, pt2])

df_seq = pd.read_csv(f"{path}/topmed_mesa_rnaseqqcv2.0.0.gene_reads.tsv.gz", sep="\t")
df_seq = df_seq.drop(columns=["Description"]).set_index("Name")

# only focus on PBMC data and visit 1
df_pt = df_pt[df_pt.visit == 5]
df_pt = df_pt[df_pt.type == "PBMC"]

# remove genes that everything is 0
df_norm = df_seq[df_pt.expr_id].transpose()
df_norm = df_norm.loc[:, ~(np.sum(df_norm, axis=0) == 0)].transpose()

# remove outliers using gtex pipline
df_rpkm = pd.read_csv(f"{path}/topmed_mesa_rnaseqqcv2.0.0.gene_rpkm.tsv.gz", sep="\t")
df_rpkm = df_rpkm.drop(columns=["Description"]).set_index("Name")
df_rpkm = df_rpkm[df_norm.columns]
df_rpkm = df_rpkm.loc[df_rpkm.index.isin(df_norm.index),:]

sample_frac_threshold = 0.2
count_threshold = 6
tpm_threshold = 0.1

df_rpkm = df_rpkm / df_rpkm.sum(0) * 1e6
ns = df_norm.shape[1]
mask = ((np.sum(df_rpkm >= tpm_threshold, axis=1) >= sample_frac_threshold * ns) & (
            np.sum(df_norm >= count_threshold, axis=1) >= sample_frac_threshold * ns)).values
df_norm = df_norm.iloc[mask, :]

df_ref = pd.read_csv("/project/nmancuso_8/data/GENCODE/gencode.v34.gene.only.tsv.gz", sep="\t")

df_norm.reset_index(inplace=True)
df_norm.Name = df_norm.Name.str.replace("\\..+", "")
df_norm = df_norm[df_norm.Name.isin(df_ref.ID2)]
# 22061 genes

# grab all the covariates
df_pt = df_pt[["expr_id", "assay_lab", "age"]]
df_pt = pd.get_dummies(df_pt, columns=["assay_lab"], drop_first=True)
df_pt.columns = ["expr_id", "age", "assay_lab"]

# read sex
pt_sex = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/freeze9/QC/sub/MESA_phs001416_TOPMed_WGS_genetic_sex.txt", sep="\t",
    names=["geno_id", "sex", "none1", "none2"], skiprows=1)

pt_sex = pt_sex[["geno_id", "sex"]]
pt_sex = pd.get_dummies(pt_sex, columns=["sex"], drop_first=True)
pt_sex.columns = ["geno_id", "sex"]

# read in ancestry
pt_ans = pd.read_csv("/project/nmancuso_8/data/MESA/processed/covariates/MESA.covar.tsv.gz", sep="\t")
pt_ans = pt_ans[["sidno", "race"]]
pt_ans.columns = ["dbgap_id", "race"]

# read in geno id
map_geno = pd.read_csv("/project/nmancuso_8/data/MESA/dbgap/WGS/freeze9/sample/sample-info.csv",
                       skiprows=10, header=None)
map_geno = map_geno[[0, 1]].rename(columns={0: "dbgap_id", 1: "geno_id"})

# read in expr id map
pt_ge_path = "/project/nmancuso_8/data/MESA/dbgap/WGS/rna_seq/sample/sample-info.csv"
pt_ge = pd.read_csv(pt_ge_path, skiprows=9, header=None)
pt_ge = pt_ge[[0, 1]]
pt_ge.columns = ["dbgap_id", "expr_id"]

df_covar = df_pt.merge(pt_ge, how="left", on="expr_id").merge(map_geno, how="inner", on="dbgap_id")\
    .merge(pt_sex, how="left", on="geno_id").merge(pt_ans, how="inner", on="dbgap_id")

dup_ids = df_covar[df_covar.dbgap_id.duplicated()].dbgap_id.unique()

rm_id_list = []
for dup_id in dup_ids:
    tmp_rm = df_covar[df_covar.dbgap_id == dup_id]
    tmp_seq = df_norm[tmp_rm.expr_id]
    csum = np.sum(tmp_seq, axis=0)
    max_idx = np.where(csum == np.max(csum))[0]
    rm_ids = csum.drop(csum[max_idx].index).index
    for rm_id in rm_ids:
        rm_id_list.append(rm_id)

df_covar = df_covar[~df_covar.expr_id.isin(rm_id_list)]

# df_covar.expr_id.duplicated().any()
# df_covar.dbgap_id.duplicated().any()
# df_covar.geno_id.duplicated().any()

# read in genotype pcs
geno_pcs = pd.read_csv(
    "/project/nmancuso_8/data/MESA/processed/covariates/genotype_pcs/all/mesa_genotype_all_pcs.eigenvec", sep="\t")
geno_pcs = geno_pcs.iloc[:, 0:11]
geno_pcs.columns = ["geno_id"] + [f"geno_pc{i}" for i in range(1, 11)]

# read in expr pcs
expr_pcs = pd.read_csv("/project/nmancuso_8/data/MESA/processed/covariates/rnaseq_pcs/rnaseq_pcs_PBMC.tsv.gz", sep="\t")
expr_pcs = expr_pcs.iloc[:, 0:16]
expr_pcs.columns = ["expr_id"] + [f"expr_pc{i}" for i in range(1, 16)]

df_all = df_covar.merge(geno_pcs, how="left", on="geno_id").merge(expr_pcs, how="left", on="expr_id")

# by ancestry
for ans in ["EUR", "AFR", "HIS"]:
    tmp_covar = df_all[df_all.race == ans]
    out_covar = tmp_covar[["geno_id", "sex", "age", "assay_lab"] + [f"geno_pc{i}" for i in range(1, 11)] + [f"expr_pc{i}" for i in range(1, 16)]]
    out_covar.to_csv(f"/project/nmancuso_8/data/MESA/processed/covariates/rnaseq/mesa_rnaseq_covar_v5_{ans}.tsv.gz", sep="\t", index=False)
    # save one copy ot scretch (to avoid bugs when running hpc that may happen from using data in the project folder)
    out_covar.to_csv(f"{scretch_path}/mesa_rnaseq_covar_v5_{ans}.tsv.gz", sep="\t", index=False, header=False)
    # output geno id for subsetting the genotype data in the analysis
    out_covar[["geno_id"]].to_csv(f"{scretch_path}/mesa_rnaseq_pt_v5_{ans}.tsv", sep="\t", index=False, header=False)
    tmp_pt = tmp_covar[["geno_id", "expr_id"]].merge(df_norm.set_index("Name").transpose().reset_index(names="expr_id"),
                                                     how="left", on="expr_id")
    tmp_pt.isna().any().any()
    tmp_pt = tmp_pt.drop(columns=["expr_id"]).set_index("geno_id")
    # this function expect the input matrix where rows are proteins and columns are participants
    tmm_counts = qtl.norm.edger_cpm(tmp_pt.transpose(), normalized_lib_sizes=True)
    tmp_pt_std = qtl.norm.inverse_normal_transform(tmm_counts).transpose().reset_index()
    tmp_pt_std.to_csv(f"/project/nmancuso_8/data/MESA/processed/expression/mesa_rnaseq_v5_{ans}.tsv.gz", sep="\t", index=False)
    tmp_pt_std.to_csv(f"{scretch_path}/mesa_rnaseq_v5_{ans}.tsv.gz", sep="\t", index=False)

tmp_pt = df_all[["geno_id", "expr_id"]].merge(df_norm.set_index("Name").transpose().reset_index(names="expr_id"),
                                                     how="left", on="expr_id")
tmp_pt.isna().any().any()
tmp_pt = tmp_pt.drop(columns=["expr_id"]).set_index("geno_id")
# this function expect the input matrix where rows are proteins and columns are participants
tmm_counts = qtl.norm.edger_cpm(tmp_pt.transpose(), normalized_lib_sizes=True)
tmp_pt_std = qtl.norm.inverse_normal_transform(tmm_counts).transpose().reset_index()
tmp_pt_std.to_csv(f"/project/nmancuso_8/data/MESA/processed/expression/mesa_rnaseq_v5_all.tsv.gz", sep="\t", index=False)
tmp_pt_std.to_csv(f"{scretch_path}/mesa_rnaseq_v5_all.tsv.gz", sep="\t", index=False)

gene_list = df_ref[df_ref.ID2.isin(df_norm.Name)].copy()
gene_list["CODE"] = df_ref.ID2 + "_" + df_ref.NAME
gene_list = gene_list.sort_values(["ID2"])

gene_list.to_csv(f"/project/nmancuso_8/data/MESA/processed/expression/mesa_rnaseq_v5_gene_list.tsv", index=False, sep="\t")
gene_list.to_csv(f"{scretch_path}/mesa_rnaseq_v5_gene_list.tsv", index=False, header=None, sep="\t")

gene_list_noMHC = gene_list[gene_list.MHC == 0]
gene_list_noMHC.to_csv(f"/project/nmancuso_8/data/MESA/processed/expression/mesa_rnaseq_v5_gene_list_noMHC.tsv", index=False, sep="\t")
gene_list_noMHC.to_csv(f"{scretch_path}/mesa_rnaseq_v5_gene_list_noMHC.tsv", index=False, header=None, sep="\t")
