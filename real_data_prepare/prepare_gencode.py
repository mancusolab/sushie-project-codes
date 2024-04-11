#!/usr/bin/env python3

import numpy as np
import pandas as pd
import qtl.io

path = "/project/nmancuso_8/data/GENCODE"
df_ref = qtl.io.gtf_to_tss_bed(f"{path}/gencode.v34.annotation.gtf.gz", feature="gene", phenotype_id="gene_id")
df_ref.reset_index(inplace=True)

# filter on chromsomes
chrs = [f'chr{i}' for i in range(1, 23)]
df_ref = df_ref.iloc[df_ref.chr.isin(chrs).values, :]

# make ensemble id
df_ref["gene_id2"] = df_ref["gene_id"].str.replace("\\..+", "")

# label MHC region chr6:28510120-33480577
df_chr6 = df_ref[df_ref.chr == "chr6"]
df_mhc = df_chr6[np.logical_and(df_chr6.end.values >= 28510120, df_chr6.start.values <= 33480577)]
df_ref["MHC"] = np.where(df_ref.index.isin(df_mhc.index), 1, 0)

chr_id = df_ref.chr.str.replace("chr", "")
df_ref["newchr"] = chr_id
df_ref["P0"] = np.where(df_ref.strand == "+", df_ref.start, df_ref.end).astype(int)
df_ref["P1"] = np.where(df_ref.strand == "+", df_ref.end, df_ref.start).astype(int)
df_ref["P0_FLANK"] = np.where(df_ref.P0 - 5e5 < 0, 0, df_ref.P0 - 5e5).astype(int)
df_ref["P1_FLANK"] = (df_ref.P1 + 5e5).astype(int)

df_ref = df_ref[["index", "gene_id2", "gene_name", "newchr", "strand", "P0", "P1", "P0_FLANK", "P1_FLANK", "MHC", "gene_type"]]
df_ref.columns = ["ID", "ID2", "NAME", "CHR", "STRAND", "TSS", "TES", "TSS_FLANK", "TES_FLANK", "MHC", "TYPE"]
df_ref.to_csv(f"{path}/gencode.v34.gene.only.tsv.gz", sep="\t", index=False)
