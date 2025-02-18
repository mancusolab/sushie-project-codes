library(tidyverse)
library(glue)

df_ref <- read_tsv(
  "/project/nmancuso_8/data/sushie/meta_data/mesa_proteins_cs_trait_list.tsv", col_names = FALSE)
lfs <- list.files("/scratch1/zeyunlu/sushie_proteins/r2")

lfs2 <- gsub("\\.cvr2\\.tsv", "", lfs)

df_out <- df_ref %>%
  filter(!X15 %in% lfs2)

write_tsv(df_out, "/scratch1/zeyunlu/fish/proteins_fish.tsv", col_names = FALSE)

df_ref <- read_tsv(
  "/project/nmancuso_8/data/sushie/meta_data/genoa_cs_trait_list.tsv", col_names = FALSE)
lfs <- list.files("/scratch1/zeyunlu/sushie_genoa/r2")

lfs2 <- gsub("\\.cvr2\\.tsv", "", lfs)

df_out <- df_ref %>%
  filter(!X2 %in% lfs2)

write_tsv(df_out, "/scratch1/zeyunlu/fish/genoa_fish.tsv", col_names = FALSE)

df_ref <- read_tsv(
  "/project/nmancuso_8/data/sushie/meta_data/mesa_rnaseq_cs_trait_list.tsv", col_names = FALSE)
lfs <- list.files("/scratch1/zeyunlu/sushie_rnaseq/r2")

lfs2 <- gsub("\\.cvr2\\.tsv", "", lfs)

df_out <- df_ref %>%
  filter(!X12 %in% lfs2)

write_tsv(df_out, "/scratch1/zeyunlu/fish/rnaseq_fish.tsv", col_names = FALSE)

