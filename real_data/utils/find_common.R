library(tidyverse)
library(glue)

scrath_paths <- "/project/nmancuso_8/data/sushie/meta_data"
df_genoa <- read_tsv(glue("{scrath_paths}/genoa_sushie_gene_list_noMHC.tsv"),
  col_names = FALSE)
df_geuv <- read_tsv(glue("{scrath_paths}/geuvadis_gene_list_noMHC.tsv"),
  col_names = FALSE)

df_same <- df_genoa %>%
  select(ID = X2) %>%
  inner_join(df_geuv %>%
      select(ID=X8),
    by = "ID")

write_tsv(df_same,
  glue("{scrath_paths}/geuvadis_overlap_gene_list_noMHC.tsv"), col_names = FALSE)

df_proteins <- read_tsv(glue("{scrath_paths}/mesa_proteins_gene_list_noMHC.tsv"),
  col_names = FALSE)
df_interval <- read_tsv(glue("{scrath_paths}/interval_gene_list_noMHC.tsv"),
  col_names = FALSE)

df_pro <- df_interval %>%
  select(ID=X8, INTERVAL=X13) %>%
  inner_join(df_proteins %>%
      select(ID=X2, MESA=X15),
    by = "ID", multiple = "all")

write_tsv(df_pro, glue("{scrath_paths}/interval_overlap_gene_list_noMHC.tsv"),
  col_names = FALSE)

df_rnaseq <- read_tsv(glue("{scrath_paths}/mesa_rnaseq_gene_list_noMHC.tsv"),
  col_names = FALSE)

df_v5 <- read_tsv(glue("{scrath_paths}/mesa_rnaseq_v5_gene_list_noMHC.tsv"),
  col_names = FALSE)

valid_v5 <- df_v5 %>%
  select(ID = X12) %>%
  inner_join(df_rnaseq %>%
      select(ID=X12))

write_tsv(valid_v5, glue("{scrath_paths}/v5_overlap_gene_list_noMHC.tsv"),
  col_names = FALSE)


