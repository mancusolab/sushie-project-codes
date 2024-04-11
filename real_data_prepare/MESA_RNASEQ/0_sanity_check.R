library(tidyverse)

setwd("/project/nmancuso_8/data/MESA/dbgap/WGS/rna_seq_v2")

df_v1 <- read_tsv("c1/TOPMed_MESA_RNAseq_Pilot_RNASeQCv2.0.0.gene_reads.gct.gz")
df_v2 <- read_tsv("c2/TOPMed_MESA_RNAseq_Pilot_RNASeQCv2.0.0.gene_reads.gct.gz")

df_total <- df_v1[,1:2774] %>%
  inner_join(df_v2[,1:173],
    by = c("Name", "Description"))

write_tsv(df_total, "topmed_mesa_rnaseqqcv2.0.0.gene_reads.tsv.gz")


df_rpkm_v1 <- read_tsv("c1/TOPMed_MESA_RNAseq_Pilot_RNASeQCv2.0.0.gene_rpkm.gct.gz")
df_rpkm_v2 <- read_tsv("c2/TOPMed_MESA_RNAseq_Pilot_RNASeQCv2.0.0.gene_rpkm.gct.gz")

df_rpkm_total <- df_rpkm_v1[,1:2774] %>%
  inner_join(df_rpkm_v2[,1:173],
    by = c("Name", "Description"))


write_tsv(df_rpkm_total, "topmed_mesa_rnaseqqcv2.0.0.gene_rpkm.tsv.gz")


