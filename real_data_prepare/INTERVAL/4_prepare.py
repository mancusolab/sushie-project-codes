#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

path = "/project/nmancuso_8/data/INTERVAL"
scratch_path = "/project/nmancuso_8/data/sushie/meta_data"

df_raw = pd.read_csv(f"{path}/phenotype/INTERVAL.proteins.plink1.tsv", sep="\t")
df_raw2 = df_raw.drop(columns=["FID"]).set_index("IID")

# 23722 proteins 462 pt
df_pc = df_raw2
df_pc -= np.mean(df_pc, axis=0)
df_pc /= np.std(df_pc, axis=0)

# calculate pcs
n_comp = 30
pca = PCA(n_components=n_comp)
raw_pcs = pca.fit_transform(df_pc)
raw_pcs -= np.mean(raw_pcs, axis=0)
raw_pcs /= np.std(raw_pcs, axis=0)
protein_pcs = pd.DataFrame(raw_pcs, index=df_pc.index, columns=[f"protein_pc{i}" for i in range(1, 31)])
protein_pcs = protein_pcs.reset_index(names="IID")
protein_pcs.to_csv(f"{path}/covariates/interval_protein_30pcs.tsv.gz", sep="\t", index=False)

# read in covariates
df_covar = pd.read_csv(f"{path}/covariates/INTERVAL.covariates.plink1.tsv", sep="\t")
df_covar = pd.get_dummies(df_covar.drop(columns="FID"), columns=["gender", "duration_blood_processing", "subchort"],
                          drop_first=True)
df_covar.columns = ["IID", "age"] + [f"geno_pcs{idx}" for idx in range(1, 4)] + ["gender", "duration_blood_processing", "subchort"]
df_covar = df_covar.merge(protein_pcs.iloc[:, 0:6], how="left", on="IID")
df_covar.to_csv(f"{path}/covariates/interval_covar.tsv.gz", sep="\t", index=False)
df_covar.to_csv(f"{scratch_path}/interval_covar.tsv", sep="\t", header=None, index=False)
df_pt = df_covar[["IID"]]
df_pt["FID"] = 0
df_pt = df_pt[["FID", "IID"]]
df_pt.to_csv(f"{scratch_path}/interval_pt_geno.tsv", sep="\t", header=None, index=False)
df_pt[["IID"]].to_csv(f"{scratch_path}/interval_pt_pheno.tsv", sep="\t", header=None, index=False)

# prepare ref file
df_ref = pd.read_csv("/project/nmancuso_8/data/GENCODE/gencode.v26lift37.GRCh37.genes.only.tsv", sep="\t")
df_ref["CHR"] = pd.to_numeric(df_ref["chr"].str.replace("chr", ""), errors="coerce")
df_ref = df_ref[~df_ref.CHR.isna()]
df_ref["CHR"] = df_ref.CHR.astype(int)
df_ref = df_ref.rename(columns={"encode_id": "encode_id_ref"})

# read in somamer id
df_soma = pd.read_csv(f"{path}/phenotype/SOMALOGIC_PROTEINS_info.tsv", sep="\t").drop(columns="TargetFullName")
df_soma = df_soma[~df_soma.UniProt.isna()]
df_soma["UniProt"] = df_soma["UniProt"].str.split(",")
df_soma = df_soma.explode("UniProt", ignore_index=True)
df_soma[["UniProt"]].to_csv("~/trash/uniprot.tsv", sep="\t", index=False, header=None)

df_map = pd.read_csv(f"{path}/phenotype/uniprot_ensemble_map.tsv", sep="\t")
df_map.columns = ["UniProt", "encode_id"]
df_map.encode_id = df_map.encode_id.str.replace("\\..+", "")
df_map = df_map.drop_duplicates()
df_soma = df_soma.merge(df_map, how="left", on="UniProt")
df_soma = df_soma.rename(columns={"encode_id": "encode_id_uni"})
# filter focus on encode
df_soma2 = df_soma.merge(df_ref, how="left", left_on="encode_id_uni", right_on="encode_id_ref")
df_soma2 = df_soma2[~df_soma2.encode_id_ref.isna()]
df_soma2 = df_soma2[["CHR", "strand", "tss", "tes", "encode_id_ref", "Target", "symbol", "UniProt", "SOMAMER_ID"]]

# focus on target
df_soma3 = df_soma.merge(df_ref, how="left", left_on="Target", right_on="symbol")
df_soma3 = df_soma3[~df_soma3.symbol.isna()]
df_soma3 = df_soma3[["CHR", "strand", "tss", "tes", "encode_id_ref", "Target", "symbol", "UniProt", "SOMAMER_ID"]]

df_final = pd.concat([df_soma2, df_soma3]).drop_duplicates()
df_final.CHR = df_final.CHR.astype(int)
df_final.tss = df_final.tss.astype(int)
df_final.tes = df_final.tes.astype(int)
df_final["TSS_FLANK"] = np.maximum(df_final.tss - 5e5, 0).astype(int)
df_final["TES_FLANK"] = (df_final.tes + 5e5).astype(int)
df_final["CODE"] = df_final.encode_id_ref + "_" + df_final.UniProt + "_" + df_final.SOMAMER_ID

# MHC region
# label MHC region GRCH37: chr6:28477797-33448354
df_chr6 = df_final[df_final.CHR == 6]
df_mhc = df_chr6[np.logical_and(df_chr6.tes.values >= 28477797, df_chr6.tss.values <= 33448354)]
df_final["MHC"] = np.where(df_final.index.isin(df_mhc.index), 1, 0)

df_final = df_final[["CHR", "strand", "tss", "tes", "TSS_FLANK", "TES_FLANK", "MHC",
                     "encode_id_ref", "Target", "symbol", "UniProt", "SOMAMER_ID", "CODE"]]
df_final.columns = ["CHR", "STRAND", "TSS", "TES", "TSS_FLANK", "TES_FLANK", "MHC",
                     "ENSEMBL_ID", "INTERVAL_SYMBOL", "GENCODE_SYMBOL", "INTERVAL_UNIPROT", "SOMAMER_ID", "CODE"]

df_final = df_final.sort_values(["CODE"])
df_final["INTERVAL_SYMBOL"] = df_final["INTERVAL_SYMBOL"].str.replace(" ", "_")

df_final.to_csv(f"{path}/interval_gene_list.tsv", index=False, sep="\t")
df_final.to_csv(f"{scratch_path}/interval_gene_list.tsv", index=False, header=None, sep="\t")

df_final_noMHC = df_final[df_final.MHC == 0]
df_final_noMHC.to_csv(f"{path}/interval_gene_list_noMHC.tsv", index=False, sep="\t")
df_final_noMHC.to_csv(f"{scratch_path}/interval_gene_list_noMHC.tsv", index=False, header=None, sep="\t")

df_raw2 = df_raw2.iloc[:, df_raw2.columns.isin(df_final.SOMAMER_ID)].reset_index(names="geno_id")
df_raw2.to_csv(f"{path}/interval_protein_levels.tsv.gz", sep="\t", index=False)
df_raw2.to_csv(f"{scratch_path}/interval_protein_levels.tsv", sep="\t", index=False)

