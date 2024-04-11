library(tidyverse)
library(broom)
library(ggpubr)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)

df_ref <- read_tsv("~/Documents/github/data/sushie_results/gencode.v34.gene.only.tsv.gz")

# EDS
df_eds <- read_excel("~/Documents/github/data/sushie_results/eds.xlsx") %>%
  dplyr::select(trait = GeneSymbol, EDS) %>%
  group_by(trait) %>%
  summarize(EDS = mean(EDS))

write_tsv(df_eds, "~/Documents/github/data/sushie_results/df_eds.tsv")

tmp_df_rvis <- read_excel("~/Documents/github/data/sushie_results/rvis.xlsx")
colnames(tmp_df_rvis) <- c("trait", "RVIS", "RVIS_percentile")

na_rvis <- tmp_df_rvis %>%
  left_join(df_ref %>%
    dplyr::select(ID2, NAME),
    by = c("trait" = "NAME")) %>%
  filter(is.na(ID2)) %>%
  dplyr::select(trait)

mapper <- bitr(na_rvis$trait, fromType="ALIAS", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

na_rvis2 <- tmp_df_rvis %>%
  left_join(df_ref %>%
      dplyr::select(ID2, NAME),
    by = c("trait" = "NAME")) %>%
  left_join(tibble(mapper),
      by = c("trait" = "ALIAS")) %>%
  mutate(trait2 = ifelse(is.na(ID2), ENSEMBL, ID2)) %>%
  filter(is.na(trait2))

traitid <- c("C10orf91", "C5orf55", "C7orf29", "CDR1", "CRIPAK",
  "FLJ42280", "FLJ43860", "GLRA4", "KIAA0664", "KIAA1045", "KIAA1804", 
  "LOC81691", "PGBD3", "PRSS45", "SPHAR", "ZHX1-C8ORF76")

ensembleid <- c("ENSG00000180066", "ENSG00000221990","ENSG00000188707",
  "ENSG00000288642", "ENSG00000179979", "ENSG00000127922",
  "ENSG00000226807", "ENSG00000188828", "ENSG00000132361",
  "ENSG00000122733", "ENSG00000143674", "ENSG00000005189", NA,
  "ENSG00000188086", "ENSG00000213029", "ENSG00000259305")

mapper2 <- tibble(traitid = traitid, eid = ensembleid)

df_rvis <- tmp_df_rvis %>%
  left_join(df_ref %>%
      dplyr::select(ID2, NAME),
    by = c("trait" = "NAME")) %>%
  left_join(tibble(mapper),
    by = c("trait" = "ALIAS")) %>%
  left_join(mapper2,
    by = c("trait" = "traitid")) %>%
  mutate(trait2 = ifelse(is.na(ID2), ENSEMBL, ID2),
    trait3 = ifelse(is.na(trait2), eid, trait2)) %>%
  filter(!is.na(trait3)) %>%
  dplyr::select(trait = trait3, RVIS) %>%
  group_by(trait) %>%
  summarize(RVIS = mean(RVIS))

write_tsv(df_rvis, "~/Documents/github/data/sushie_results/df_rvis.tsv")

tmp_df_pli <- read_excel("~/Documents/github/data/sushie_results/pli.xlsx", sheet = 2) %>%
  dplyr::select(gene, pLI)

mapper <- clusterProfiler::bitr(tmp_df_pli$gene, fromType="ALIAS", toType="ENSEMBL", OrgDb="org.Hs.eg.db")

df_pli1 <- tmp_df_pli %>%
  left_join(tibble(mapper),
    by = c("gene" = "ALIAS")) %>%
  filter(!is.na(ENSEMBL)) %>%
  dplyr::select(trait = ENSEMBL, pLI)

na_pli <- tmp_df_pli %>%
  left_join(tibble(mapper),
    by = c("gene" = "ALIAS")) %>%
  filter(is.na(ENSEMBL))

df_pli2 <- na_pli %>%
  mutate(gene = tolower(gene)) %>%
  left_join(df_ref %>%
      mutate(NAME = tolower(NAME)) %>%
      dplyr::select(NAME, ID2),
    by = c("gene" = "NAME")) %>%
  filter(!is.na(ID2)) %>%
  dplyr::select(trait = ID2, pLI)
  

df_pli <- bind_rows(df_pli1, df_pli2) %>%
  group_by(trait) %>%
  summarize(pLI = mean(pLI))

write_tsv(df_pli, "~/Documents/github/data/sushie_results/df_pli.tsv")

df_shet <- read_tsv("~/Documents/github/data/sushie_results/s_het.tsv") %>%
  dplyr::select(trait = ensg, s_het=post_mean)
  
write_tsv(df_shet, "~/Documents/github/data/sushie_results/df_shet.tsv")
  

df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/gencode.v34.gene.only.tsv.gz")

df_tmp <- read_tsv("~/Downloads/gnomad.v4.0.constraint_metrics.tsv") %>%
  filter(grepl("ENST", transcript)) %>%
  select(gene, pLI = lof.pLI, LOEUF = lof.oe_ci.upper, Z = syn.z_score) %>%
  pivot_longer(cols = c(pLI, LOEUF)) %>%
  group_by(gene, name) %>%
  summarize(value= mean(value, na.rm = TRUE))

new_pli <- df_tmp %>%
  inner_join(df_ref, by = c("gene" = "NAME")) %>%
  ungroup() %>%
  select(ID2, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  select(trait = ID2, LOEUF, pLI) %>%
  pivot_longer(cols=c(2:3)) %>%
  filter(!is.na(value))

write_tsv(new_pli, "~/Documents/github/data/sushie_results/Constraint/df_pli_new.tsv")

