library(tidyverse)
library(broom)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
real_data_path <- "~/Documents/github/data/sushie_results/real2/"
meta_data_path <- "your_path_from_zenodo"

df_snps <- read_tsv(glue("{real_data_path}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  filter(pip_all > 0.9) %>%
  select(trait, snp) %>%
  mutate(method = "SuShiE") %>%
  bind_rows(
    read_tsv(glue("{real_data_path}/rnaseq_susiex_cs.tsv.gz")) %>%
      filter(!is.na(snp)) %>%
      filter(pip_all > 0.9) %>%
      mutate(method = "SuSiEx") %>%
      select(trait, snp, method),
    read_tsv(glue("{real_data_path}/rnaseq_mesusie_cs.tsv.gz")) %>%
      filter(!is.na(snp)) %>%
      filter(pip_all > 0.9) %>%
      mutate(method = "MESuSiE") %>%
      select(trait, snp, method)
  ) %>%
  group_by(trait, snp) %>%
  filter(n()==3) %>%
  ungroup() %>%
  distinct(trait, snp)

df_snps %>%
  distinct(trait)

df_snps %>%
  distinct(snp)

df_meta <- read_tsv(glue("{meta_data_path}/mesa_rnaseq_gene_list_noMHC.tsv"),
  col_names = FALSE) %>%
  inner_join(df_snps %>%
      distinct(trait),
    by = c("X12"="trait"))

# write_tsv(df_meta, "~/Documents/github/sushie-data-codes/sim_codes/param/real_qtl.tsv", col_names = FALSE)

for (gene_name in unique(df_meta$X12)){
  df_tmp <- df_snps %>%
    filter(trait == gene_name) %>%
    select(snp)
  # write_tsv(df_tmp, paste0("~/Documents/github/sushie-data-codes/sim_codes/param/gene_level/real_qtl_", gene_name, ".tsv"), col_names = FALSE)
}


