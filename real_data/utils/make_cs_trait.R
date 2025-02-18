library(tidyverse)
library(glue)
library(broom)

trait_list <- read_tsv("~/data/sushie/real2/rnaseq_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  bind_rows(
    read_tsv("~/data/sushie/real2/rnaseq_indep.sushie_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/rnaseq_normal.meta_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/rnaseq_normal.mega_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/rnaseq_susiex_cs.tsv.gz") %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/rnaseq_mesusie_cs.tsv.gz") %>%
      distinct(trait)
  ) %>%
  distinct(trait)

df_ref <- read_tsv("/project/nmancuso_8/data/sushie/meta_data/mesa_rnaseq_gene_list_noMHC.tsv",
  col_names=FALSE)
new_df_ref <- df_ref %>%
  filter(X12 %in% trait_list$trait)

write_tsv(new_df_ref, "/project/nmancuso_8/data/sushie/meta_data/mesa_rnaseq_cs_trait_list.tsv", col_names = FALSE)


trait_list <- read_tsv("~/data/sushie/real2/proteins_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  bind_rows(
    read_tsv("~/data/sushie/real2/proteins_indep.sushie_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/proteins_normal.meta_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/proteins_normal.mega_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/proteins_susiex_cs.tsv.gz") %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/proteins_mesusie_cs.tsv.gz") %>%
      distinct(trait)
  ) %>%
  distinct(trait)

df_ref <- read_tsv("/project/nmancuso_8/data/sushie/meta_data/mesa_proteins_gene_list_noMHC.tsv",
  col_names=FALSE)
new_df_ref <- df_ref %>%
  filter(X15 %in% trait_list$trait)

write_tsv(new_df_ref, "/project/nmancuso_8/data/sushie/meta_data/mesa_proteins_cs_trait_list.tsv", col_names = FALSE)

trait_list <- read_tsv("~/data/sushie/real2/genoa_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  bind_rows(
    read_tsv("~/data/sushie/real2/genoa_indep.sushie_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/genoa_normal.meta_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/genoa_normal.mega_cs.tsv.gz") %>%
      filter(!is.na(snp)) %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/genoa_susiex_cs.tsv.gz") %>%
      distinct(trait),
    
    read_tsv("~/data/sushie/real2/genoa_mesusie_cs.tsv.gz") %>%
      distinct(trait)
  ) %>%
  distinct(trait)

df_ref <- read_tsv("/project/nmancuso_8/data/sushie/meta_data/genoa_sushie_gene_list_noMHC.tsv",col_names=FALSE)
new_df_ref <- df_ref %>%
  filter(X2 %in% trait_list$trait)

write_tsv(new_df_ref, "/project/nmancuso_8/data/sushie/meta_data/genoa_cs_trait_list.tsv", col_names = FALSE)
