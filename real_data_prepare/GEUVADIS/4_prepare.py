#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

path = "/project/nmancuso_8/data/GEUVADIS/ProcessedData"
scretch_path = "/project/nmancuso_8/data/sushie/meta_data"

df_raw = pd.read_csv(f"{path}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", sep="\t")
df_raw2 = df_raw.drop(columns=["Gene_Symbol", "Coord", "Chr"]).set_index("TargetID")

# 23722 proteins 462 pt
df_pc = df_raw2.transpose()
df_pc -= np.mean(df_pc, axis=0)
df_pc /= np.std(df_pc, axis=0)

# calculate pcs
n_comp = 30
pca = PCA(n_components=n_comp)
raw_pcs = pca.fit_transform(df_pc)
raw_pcs -= np.mean(raw_pcs, axis=0)
raw_pcs /= np.std(raw_pcs, axis=0)
expr_pcs = pd.DataFrame(raw_pcs, index=df_pc.index, columns=[f"expr_pc{i}" for i in range(1, 31)])
expr_pcs = expr_pcs.reset_index(names="IID")
expr_pcs.to_csv(f"{path}/geuvadis_expr_30pcs.tsv.gz", sep="\t", index=False)

# remove low expressed genes
sample_frac_threshold = 0.2
tpm_threshold = 0.1
df_rpkm = df_raw2
df_rpkm = df_rpkm / df_rpkm.sum(0) * 1e6
ns = df_rpkm.shape[1]
mask = (np.sum(df_rpkm >= tpm_threshold, axis=1) >= sample_frac_threshold * ns).values
df_raw2 = df_raw2.iloc[mask, :]

df_raw2.index = df_raw2.index.str.replace("\\..+", "")

# prepare reference file
df_ref = pd.read_csv("/project/nmancuso_8/data/GENCODE/gencode.v26lift37.GRCh37.genes.only.tsv", sep="\t")
df_ref["CHR"] = pd.to_numeric(df_ref["chr"].str.replace("chr", ""), errors="coerce")
df_ref = df_ref[~df_ref.CHR.isna()]
df_ref["CHR"] = df_ref.CHR.astype(int)

df_ref = df_ref[df_ref.encode_id.isin(df_raw2.index)]
df_raw2 = df_raw2[df_raw2.index.isin(df_ref.encode_id)]

# MHC region
# label MHC region GRCH37: chr6:28477797-33448354
df_chr6 = df_ref[df_ref.chr == "chr6"]
df_mhc = df_chr6[np.logical_and(df_chr6.tes.values >= 28477797, df_chr6.tss.values <= 33448354)]
df_ref["MHC"] = np.where(df_ref.index.isin(df_mhc.index), 1, 0)

df_ref["TSS_FLANK"] = np.maximum(df_ref.tss - 5e5, 0).astype(int)
df_ref["TES_FLANK"] = (df_ref.tes + 5e5).astype(int)

df_ref = df_ref[["CHR", "strand", "tss", "tes", "TSS_FLANK", "TES_FLANK", "MHC", "encode_id", "symbol"]]
df_ref.columns = ["CHR", "STRAND", "TSS", "TES",  "TSS_FLANK", "TES_FLANK", "MHC", "ID", "NAME"]

df_ref = df_ref.sort_values(["ID"])

df_ref.to_csv(f"{path}/geuvadis_gene_list.tsv", index=False, sep="\t")
df_ref.to_csv(f"{scretch_path}/geuvadis_gene_list.tsv", index=False, header=None, sep="\t")

df_ref_noMHC = df_ref[df_ref.MHC == 0]
df_ref_noMHC.to_csv(f"{path}/geuvadis_gene_list_noMHC.tsv", index=False, sep="\t")
df_ref_noMHC.to_csv(f"{scretch_path}/geuvadis_gene_list_noMHC.tsv", index=False, header=None, sep="\t")

# prepare covariate
# read in expr pcs
expr_pcs = pd.read_csv(f"{path}/geuvadis_expr_30pcs.tsv.gz", sep="\t")

geno_pcs = pd.read_csv("/project/nmancuso_8/data/GEUVADIS/ProcessedData/geno_pcs/geuvadis_genotype_all_pcs.eigenvec", sep="\t")
geno_pcs = geno_pcs.drop(columns="#FID")
geno_pcs.columns = ["IID"] + [f"geno_pc{idx}" for idx in range(1, 31)]

for pop in ["EUR", "YRI"]:
    df_ans = pd.read_csv(f"/project/nmancuso_8/data/GEUVADIS/Populations/original/{pop}_sex.tsv", header=None, sep="\t")
    df_ans = df_ans[[1, 2]]
    df_ans.columns = ["IID", "SEX"]
    df_covar = df_ans.merge(expr_pcs.iloc[:, 0:6], how="left", on="IID").merge(geno_pcs.iloc[:, 0:6], how="left", on="IID")
    df_covar.to_csv(f"{path}/geuvadis_{pop}_covar.tsv.gz", sep="\t", index=False)
    df_covar.to_csv(f"{scretch_path}/geuvadis_{pop}_covar.tsv", sep="\t", header=None, index=False)
    df_covar[["IID", "IID"]].to_csv(f"{scretch_path}/geuvadis_{pop}_pt.tsv", sep="\t", header=None, index=False)
    df_pheno = df_raw2.transpose().reset_index(names="IID")
    df_pheno = df_pheno[df_pheno.IID.isin(df_covar.IID)]
    df_pheno = df_pheno.rename(columns={"IID": "geno_id"})
    df_pheno.to_csv(f"{path}/geuvadis_{pop}_expression_rpkm.tsv.gz", sep="\t", index=False)
    df_pheno.to_csv(f"{scretch_path}/geuvadis_{pop}_expression_rpkm.tsv", sep="\t", index=False)

