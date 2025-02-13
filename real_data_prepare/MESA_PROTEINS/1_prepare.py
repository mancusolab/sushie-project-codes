#!/usr/bin/env python3

import pandas as pd
import numpy as np
import qtl.norm
from sklearn.decomposition import PCA
import copy

path = "/project/nmancuso_8/data/MESA/dbgap/WGS/proteome/"
scretch_path = "/project/nmancuso_8/data/sushie/meta_data"

df_pt1 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010511.v1.p1.c1.TOPMed_MESA_Proteomics_Sample_Attributes.HMB.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "pro_id",
                                           "SUBJECT_ID": "dbgap_id",
                                           "AGE_AT_COLLECTION": "age",
                                            "COLLECTION_VISIT": "visit"})
df_pt1 = df_pt1[["pro_id", "dbgap_id", "age", "visit"]]

df_pt2 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010511.v1.p1.c2.TOPMed_MESA_Proteomics_Sample_Attributes.HMB-NPU.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "pro_id",
                                           "SUBJECT_ID": "dbgap_id",
                                           "AGE_AT_COLLECTION": "age",
                                           "COLLECTION_VISIT": "visit"
                                           })

df_pt2 = df_pt2[["pro_id", "dbgap_id", "age", "visit"]]
df_pt = pd.concat([df_pt1, df_pt2])
v1_pt = df_pt[df_pt.visit == 1]

df_c1 = pd.read_csv(f"{path}/c1/MESA_ProteomicsDataMergedRunlistKey_DS_HMB_20220426.txt", sep="\t")
df_c2 = pd.read_csv(f"{path}/c2/MESA_ProteomicsDataMergedRunlistKey_DS_HMB-NPU_20220426.txt", sep="\t")
df_pro = np.log(pd.concat([df_c1, df_c2]).set_index("TOP_ID").transpose())[v1_pt.pro_id.values]
df_pro_std = qtl.norm.inverse_normal_transform(df_pro).transpose()
# 1317 proteins 1966 pt

# run pcs on visit 1 measurements
# scale first
df_pro_std -= np.mean(df_pro_std, axis=0)
df_pro_std /= np.std(df_pro_std, axis=0)

n_comp = 30
pca = PCA(n_components=n_comp)
raw_pcs = pca.fit_transform(df_pro_std)
raw_pcs -= np.mean(raw_pcs, axis=0)
raw_pcs /= np.std(raw_pcs, axis=0)
proteins_pcs = pd.DataFrame(raw_pcs, index=df_pro_std.index, columns=[f"proteins_pc{i}" for i in range(1, 31)])
proteins_pcs = proteins_pcs.reset_index(names="pro_id")
proteins_pcs.to_csv("/project/nmancuso_8/data/MESA/processed/covariates/protein_pcs/mesa_proteins_all_pcs_v1.tsv.gz",
                  sep="\t", index=False)

df_pro = pd.concat([df_c1, df_c2])

# prepare protein file
df_ref = pd.read_csv("/project/nmancuso_8/data/GENCODE/gencode.v34.gene.only.tsv.gz", sep="\t")
df_ref2 = df_ref[["ID2", "NAME"]]
df_ref3 = copy.deepcopy(df_ref2)
df_ref3.columns = ["ID2", "EntrezGeneSymbol"]
df_refmap = pd.read_csv(f"/project/nmancuso_8/data/MESA/processed/protein/aptamer_map.tsv", sep="\t")
df_refmap.columns = ["ID2", "SomaId"]
pro_name = pd.read_csv(f"{path}/c1/MESA_ProteomicsDataMergedRunlistKey_ProteinInfo_All_20220426.txt", sep="\t")
pro_name = pro_name[["SomaId", "Target", "EntrezGeneSymbol"]]
tmp_ref = pro_name.merge(df_refmap, how="left", on="SomaId").merge(df_ref2, how="left", on="ID2").\
    merge(df_ref3, how="left", on="EntrezGeneSymbol")

f_tmp_ref1 = tmp_ref[tmp_ref["ID2_x"] == tmp_ref["ID2_y"]]
f_tmp_ref1["ID2"] = f_tmp_ref1["ID2_x"]
f_tmp_ref1 = f_tmp_ref1[["SomaId", "Target", "EntrezGeneSymbol", "ID2"]]

tmp_ref2 = tmp_ref[~(tmp_ref["ID2_x"] == tmp_ref["ID2_y"])]

f_tmp_ref2 = tmp_ref2[~(tmp_ref2["ID2_x"].notna() & tmp_ref2["ID2_y"].notna())]
f_tmp_ref2["ID2"] = np.where(f_tmp_ref2["ID2_x"].notna(), f_tmp_ref2["ID2_x"], f_tmp_ref2["ID2_y"])
f_tmp_ref2 = f_tmp_ref2[~f_tmp_ref2.ID2.isna()]
f_tmp_ref2 = f_tmp_ref2[["SomaId", "Target", "EntrezGeneSymbol", "ID2"]]

f_tmp_ref3 = tmp_ref2[(tmp_ref2["ID2_x"].notna() & tmp_ref2["ID2_y"].notna())]
f_tmp_ref3["ID2"] = f_tmp_ref3.ID2_x
f_tmp_ref3 = f_tmp_ref3[["SomaId", "Target", "EntrezGeneSymbol", "ID2"]]

f_tmp_ref4 = tmp_ref2[(tmp_ref2["ID2_x"].notna() & tmp_ref2["ID2_y"].notna())]
f_tmp_ref4 = f_tmp_ref4[["SomaId", "Target", "NAME", "ID2_y"]]
f_tmp_ref4.columns = ["SomaId", "Target", "EntrezGeneSymbol", "ID2"]

f_tmp_ref = pd.concat([f_tmp_ref1, f_tmp_ref2, f_tmp_ref3, f_tmp_ref4])
f_tmp_ref = f_tmp_ref.drop_duplicates()
f_tmp_ref["CODE"] = f_tmp_ref.ID2 + "_" + f_tmp_ref.SomaId

df_pt1 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010511.v1.p1.c1.TOPMed_MESA_Proteomics_Sample_Attributes.HMB.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "pro_id",
                                           "SUBJECT_ID": "dbgap_id",
                                           "AGE_AT_COLLECTION": "age",
                                            "COLLECTION_VISIT": "visit"})
df_pt1 = df_pt1[["pro_id", "dbgap_id", "age", "visit"]]

df_pt2 = pd.read_csv(
    "/project/nmancuso_8/data/MESA/dbgap/WGS/phs001416.v3.pht010511.v1.p1.c2.TOPMed_MESA_Proteomics_Sample_Attributes.HMB-NPU.txt.gz",
    sep="\t", skiprows=10).rename(columns={"SAMPLE_ID": "pro_id",
                                           "SUBJECT_ID": "dbgap_id",
                                           "AGE_AT_COLLECTION": "age",
                                           "COLLECTION_VISIT": "visit"
                                           })

df_pt2 = df_pt2[["pro_id", "dbgap_id", "age", "visit"]]
df_pt = pd.concat([df_pt1, df_pt2])
# only use first visit
df_pt = df_pt[df_pt.visit == 1]

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

df_all = pt_ans.merge(df_pt, how="inner", on="dbgap_id").merge(map_geno, how="inner", on="dbgap_id").merge(pt_sex, how="inner", on="geno_id")
# df_all.dbgap_id.duplicated().any() false, no duplicated id
# df_all.pro_id.duplicated().any() false, no duplicated id
# df_all.geno_id.duplicated().any() false, no duplicated id

# read in genotype pcs
geno_pcs = pd.read_csv(
    "/project/nmancuso_8/data/MESA/processed/covariates/genotype_pcs/all/mesa_genotype_all_pcs.eigenvec", sep="\t")
geno_pcs = geno_pcs.iloc[:, 0:11]
geno_pcs.columns = ["geno_id"] + [f"geno_pc{i}" for i in range(1, 11)]

df_covar = df_all.merge(geno_pcs, how="left", on="geno_id").merge(proteins_pcs, how="left", on="pro_id")

# all participants
out_covar = df_covar[["geno_id", "dbgap_id", "pro_id", "sex", "age"] + [f"geno_pc{i}" for i in range(1, 11)] + [f"proteins_pc{i}" for i in range(1, 16)]]
# output everyone with headers
out_covar.to_csv(f"/project/nmancuso_8/data/MESA/processed/covariates/proteins/mesa_proteins_covar_v1_all.tsv.gz", sep="\t", index=False)

tmp_pt = out_covar[["geno_id", "pro_id"]].merge(df_pro, how="left", left_on="pro_id", right_on="TOP_ID")
tmp_pt = tmp_pt.drop(columns=["pro_id", "TOP_ID"]).set_index("geno_id")
tmp_pt_std = qtl.norm.inverse_normal_transform(np.log(tmp_pt).transpose()).transpose().reset_index()
melted_df = pd.melt(tmp_pt_std, id_vars=["geno_id"], var_name="SomaId", value_name="value")
melted_df = melted_df.merge(f_tmp_ref, how="left", on="SomaId").drop(columns=["SomaId", "ID2", "Target", "EntrezGeneSymbol"])
melted_df = melted_df[~melted_df.CODE.isna()]
melted_df = melted_df.pivot(index="geno_id", columns="CODE", values="value").reset_index()
melted_df.to_csv(f"/project/nmancuso_8/data/MESA/processed/protein/mesa_proteins_v1_all.tsv.gz", sep="\t", index=False)

# by ancestry
for ans in ["EUR", "AFR", "HIS", "EAS"]:
    tmp_covar = df_covar[df_covar.race == ans]
    out_covar = tmp_covar[["geno_id", "sex", "age"] + [f"geno_pc{i}" for i in range(1, 11)] + [f"proteins_pc{i}" for i in range(1, 16)]]
    out_covar.to_csv(f"/project/nmancuso_8/data/MESA/processed/covariates/proteins/mesa_proteins_covar_v1_{ans}.tsv.gz", sep="\t", index=False)
    # save one copy ot scretch (to avoid bugs when running hpc that may happen from using data in the project folder)
    out_covar.to_csv(f"{scretch_path}/mesa_proteins_covar_{ans}.tsv.gz", sep="\t", index=False, header=False)
    # output geno id for subsetting the genotype data in the analysis
    out_covar[["geno_id"]].to_csv(f"{scretch_path}/mesa_proteins_pt_{ans}.tsv", sep="\t", index=False, header=False)
    tmp_pt = tmp_covar[["geno_id", "pro_id"]].merge(df_pro, how="left", left_on="pro_id", right_on="TOP_ID")
    tmp_pt.isna().any().any()
    tmp_pt = tmp_pt.drop(columns=["pro_id", "TOP_ID"]).set_index("geno_id")
    # this function expect the input matrix where rows are proteins and columns are participants
    tmp_pt_std = qtl.norm.inverse_normal_transform(np.log(tmp_pt).transpose()).transpose().reset_index()
    melted_df = pd.melt(tmp_pt_std, id_vars=["geno_id"], var_name="SomaId", value_name="value")
    melted_df = melted_df.merge(f_tmp_ref, how="left", on="SomaId").drop(
        columns=["SomaId", "ID2", "Target", "EntrezGeneSymbol"])
    melted_df = melted_df[~melted_df.CODE.isna()]
    melted_df = melted_df.pivot(index="geno_id", columns="CODE", values="value").reset_index()
    melted_df.to_csv(f"/project/nmancuso_8/data/MESA/processed/protein/mesa_proteins_v1_{ans}.tsv.gz", sep="\t", index=False)
    melted_df.to_csv(f"{scretch_path}/mesa_proteins_v1_{ans}.tsv.gz", sep="\t", index=False)

df_final = df_ref.merge(f_tmp_ref, how="right", on="ID2")
df_final["Target"] = df_final["Target"].str.replace(r"\s+", "_")
df_final["EntrezGeneSymbol"] = df_final["EntrezGeneSymbol"].str.replace(r"\s+", "_")
df_final = df_final.sort_values(["ID2"])

df_final.to_csv(f"/project/nmancuso_8/data/MESA/processed/protein/mesa_proteins_gene_list.tsv", index=False, sep="\t")
df_final.to_csv(f"{scretch_path}/mesa_proteins_gene_list.tsv", index=False, header=None, sep="\t")

df_final_noMHC = df_final[df_final.MHC == 0]
df_final_noMHC.to_csv(f"/project/nmancuso_8/data/MESA/processed/protein/mesa_proteins_gene_list_noMHC.tsv", index=False, sep="\t")
df_final_noMHC.to_csv(f"{scretch_path}/mesa_proteins_gene_list_noMHC.tsv", index=False, header=None, sep="\t")
# 1274 genes

# make sure no duplicate genotype
df_dup = pd.read_csv("/project/nmancuso_8/data/MESA/dbgap/WGS/freeze9/sample/sample-info.csv", skiprows=10, header=None)
(df_dup[5] == "N").all()
