library(tidyverse)
library(ggpubr)
library(broom)
source("./utils.R")

# mesa.rnaseq
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz")

rnaseq_indep <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_indep.sushie_cs.tsv.gz")

rnaseq_meta <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.meta_cs.tsv.gz")

rnaseq_mega <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.mega_cs.tsv.gz")

# v5
v5_cov <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.sushie_cs.tsv.gz")

v5_indep <- read_tsv("~/Documents/github/data/sushie_results/real/v5_indep.sushie_cs.tsv.gz")

v5_meta <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.meta_cs.tsv.gz")

v5_mega <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.mega_cs.tsv.gz")

# mesa.proteins
proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz")

proteins_indep <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_indep.sushie_cs.tsv.gz")

proteins_meta <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.meta_cs.tsv.gz")

proteins_mega <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.mega_cs.tsv.gz")

# interval
interval_cov <- read_tsv("~/Documents/github/data/sushie_results/real/interval_normal.sushie_cs.tsv.gz")


# genoa
genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz")

genoa_indep <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_indep.sushie_cs.tsv.gz")

genoa_meta <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.meta_cs.tsv.gz")

genoa_mega <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.mega_cs.tsv.gz")


# geuvadis
geuvadis_cov <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.sushie_cs.tsv.gz")

geuvadis_indep <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_indep.sushie_cs.tsv.gz")

geuvadis_meta <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.meta_cs.tsv.gz")

geuvadis_mega <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.mega_cs.tsv.gz")

# her
rnaseq_her <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_her.tsv.gz")

proteins_her <-
  read_tsv("~/Documents/github/data/sushie_results/real/proteins_her.tsv.gz")

genoa_her <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_her.tsv.gz")

v5_her <- read_tsv("~/Documents/github/data/sushie_results/real/v5_her.tsv.gz")

interval_her <-
  read_tsv("~/Documents/github/data/sushie_results/real/interval_her.tsv.gz")

geuvadis_her <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_her.tsv.gz")

# valid
rnaseq_valid_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_valid.normal.tsv.gz")

proteins_valid_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_valid.normal.tsv.gz")

genoa_valid_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_valid.normal.tsv.gz")

rnaseq_valid_indep <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_valid.indep.tsv.gz")

rnaseq_valid_meta <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_valid.meta.tsv.gz")

rnaseq_valid_mega <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_valid.mega.tsv.gz")

genoa_valid_indep <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_valid.indep.tsv.gz")

genoa_valid_meta <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_valid.meta.tsv.gz")

genoa_valid_mega <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_valid.mega.tsv.gz")


df_ref_overlap_rnaseq <- read_tsv("~/Documents/github/data/sushie_results/metadata/v5_overlap_gene_list_noMHC.tsv", col_names = FALSE)
df_ref_overlap_proteins <- read_tsv("~/Documents/github/data/sushie_results/metadata/interval_overlap_gene_list_noMHC.tsv", col_names = FALSE)
df_ref_overlap_genoa <- read_tsv("~/Documents/github/data/sushie_results/metadata/geuvadis_overlap_gene_list_noMHC.tsv", col_names = FALSE)

rnaseq_sig_trait <- rnaseq_her %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==3) %>%
  distinct(trait)

proteins_sig_trait <- proteins_her %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==3) %>%
  distinct(trait)

genoa_sig_trait <- genoa_her %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==2) %>%
  distinct(trait)

v5_sig_trait <- v5_her %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==3) %>%
  distinct(trait)

interval_sig_trait <- interval_her %>%
  filter(p_value < 0.05) %>%
  distinct(trait)

geuvadis_sig_trait <- geuvadis_her %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==2) %>%
  distinct(trait)

overlap_rnaseq_v5 <- rnaseq_cov %>%
  filter(trait %in% df_ref_overlap_rnaseq$X1) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  inner_join(v5_cov %>%
      distinct(trait))

overlap_proteins_interval <- df_ref_overlap_proteins %>%
  inner_join(proteins_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait) %>%
      select(X3 = trait)) %>%
  inner_join(interval_cov %>%
      distinct(trait) %>%
      select(X2 = trait)) %>%
  select(trait = X3)

overlap_genoa_geuvadis <- genoa_cov %>%
  filter(trait %in% df_ref_overlap_genoa$X1) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  inner_join(geuvadis_cov %>%
      distinct(trait))

nrow(overlap_rnaseq_v5) + nrow(overlap_proteins_interval) +
  nrow(overlap_genoa_geuvadis)

# across all studies
rnaseq_cov %>%
  filter(trait %in% overlap_rnaseq_v5$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac, z_cos)) %>%
  bind_rows(proteins_cov %>%
      filter(trait %in% overlap_proteins_interval$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(proteins_valid_cov %>%
          filter(match == 1) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac, z_cos)),
    genoa_cov %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_cov %>%
          filter(match == 1) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac, z_cos))) %>%
  mutate(z_cos = ifelse(is.infinite(z_cos), max(z_cos),z_cos)) %>%
  summarize(dem = n(),
    num = sum(match, na.rm = TRUE),
    ratio = num/dem,
    m_zcos = mean(z_cos, na.rm = TRUE),
    m_cos = mean(cos, na.rm = TRUE),
    m_wjac = mean(w_jac, na.rm = TRUE),
    m_nhac = mean(n_jac, na.rm = TRUE))

# rnaseq
rnaseq_cov %>%
  filter(trait %in% overlap_rnaseq_v5$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  summarize(dem = n(),
    num = sum(match, na.rm = TRUE),
    ratio = num/dem,
    m_cos = mean(cos, na.rm = TRUE),
    m_wjac = mean(w_jac, na.rm = TRUE),
    m_nhac = mean(n_jac, na.rm = TRUE))

# proteins
proteins_cov %>%
  filter(trait %in% overlap_proteins_interval$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(proteins_valid_cov %>%
      filter(match == 1) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  summarize(dem = n(),
    num = sum(match, na.rm = TRUE),
    ratio = num/dem,
    m_cos = mean(cos, na.rm = TRUE),
    m_wjac = mean(w_jac, na.rm = TRUE),
    m_nhac = mean(n_jac, na.rm = TRUE))

# genoa
genoa_cov %>%
  filter(trait %in% overlap_genoa_geuvadis$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(genoa_valid_cov %>%
      filter(match == 1) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  summarize(dem = n(),
    num = sum(match, na.rm = TRUE),
    ratio = num/dem,
    m_cos = mean(cos, na.rm = TRUE),
    m_wjac = mean(w_jac, na.rm = TRUE),
    m_nhac = mean(n_jac, na.rm = TRUE))


rnaseq_cov %>%
  filter(CSIndex == 1) %>%
  filter(trait %in% overlap_rnaseq_v5$trait) %>%
  filter(trait %in% rnaseq_sig_trait$trait) %>%
  filter(trait %in% v5_sig_trait$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(proteins_cov %>%
      filter(CSIndex == 1) %>%
      filter(trait %in% proteins_sig_trait$trait) %>%
      filter(trait %in% interval_sig_trait$trait) %>%
      filter(trait %in% overlap_proteins_interval$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(proteins_valid_cov %>%
          filter(match == 1) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)),
    genoa_cov %>%
      filter(CSIndex == 1) %>%
      filter(trait %in% genoa_sig_trait$trait) %>%
      filter(trait %in% geuvadis_sig_trait$trait) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_cov %>%
          filter(match == 1) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) %>%
  summarize(dem = n(),
    num = sum(match, na.rm = TRUE),
    ratio = num/dem,
    m_cos = mean(cos, na.rm = TRUE),
    m_wjac = mean(w_jac, na.rm = TRUE),
    m_nhac = mean(n_jac, na.rm = TRUE))

rep_comp_ratio <- tibble()
rep_comp_cos <- tibble()

# cov and indep
cov_indep_genes_rnaseq <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(rnaseq_indep %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

cov_indep_genes_genoa <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(genoa_indep %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

# sushie cov
cov_indep_res <- rnaseq_cov %>%
  inner_join(cov_indep_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait)%>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_cov %>%
      inner_join(cov_indep_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_cov %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 

# sushie indep
indep_res <- rnaseq_indep %>%
  inner_join(cov_indep_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_indep %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_indep %>%
      inner_join(cov_indep_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_indep %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 

sum(cov_indep_res$match, na.rm = TRUE)
sum(indep_res$match, na.rm = TRUE)

rep_comp_ratio <- rep_comp_ratio %>%
  bind_rows(
    tidy(prop.test(c(sum(cov_indep_res$match, na.rm = TRUE),
      sum(indep_res$match, na.rm = TRUE)),
      c(nrow(cov_indep_res), nrow(indep_res)))) %>%
      mutate(type = "Cov VS Indep")
  )

mean(cov_indep_res$cos, na.rm = TRUE)
mean(indep_res$cos, na.rm = TRUE)

rep_comp_cos <- rep_comp_cos %>%
  bind_rows(
    tidy(t.test(cov_indep_res$cos, indep_res$cos)) %>%
      mutate(type = "Cov VS Indep"))

# cov and meta
cov_meta_genes_rnaseq <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(rnaseq_meta %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

cov_meta_genes_genoa <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(genoa_meta %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

# sushie cov
cov_meta_res <- rnaseq_cov %>%
  inner_join(cov_meta_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait)%>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_cov %>%
      inner_join(cov_meta_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_cov %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 

# sushie meta
meta_res <- rnaseq_meta %>%
  inner_join(cov_meta_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait) %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_meta %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      distinct(main_gene, main_cs_index, match, cos, w_jac, n_jac) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_meta %>%
      inner_join(cov_meta_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_meta %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 

sum(cov_meta_res$match, na.rm = TRUE)
sum(meta_res$match, na.rm = TRUE)

rep_comp_ratio <- rep_comp_ratio %>%
  bind_rows(
    tidy(prop.test(c(sum(cov_meta_res$match, na.rm = TRUE),
      sum(meta_res$match, na.rm = TRUE)),
      c(nrow(cov_meta_res), nrow(meta_res)))) %>%
      mutate(type = "Cov VS Meta"))

mean(cov_meta_res$cos, na.rm = TRUE)
mean(meta_res$cos, na.rm = TRUE)
rep_comp_cos <- rep_comp_cos %>%
  bind_rows(
    tidy(t.test(cov_meta_res$cos, meta_res$cos)) %>%
      mutate(type = "Cov VS Meta"))


# cov and mega
cov_mega_genes_rnaseq <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(rnaseq_mega %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

cov_mega_genes_genoa <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  inner_join(genoa_mega %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex))

# sushie cov
cov_mega_res <- rnaseq_cov %>%
  inner_join(cov_mega_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait)%>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_cov %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_cov %>%
      inner_join(cov_mega_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_cov %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 


# sushie mega
mega_res <- rnaseq_mega %>%
  inner_join(cov_mega_genes_rnaseq) %>%
  filter(trait %in% overlap_rnaseq_v5$trait)%>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  left_join(rnaseq_valid_mega %>%
      filter(match == 1) %>%
      group_by(main_gene, main_cs_index) %>%
      filter(cos == max(cos)) %>%
      select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac)) %>%
  bind_rows(genoa_mega %>%
      inner_join(cov_mega_genes_genoa) %>%
      filter(trait %in% overlap_genoa_geuvadis$trait) %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex) %>%
      left_join(genoa_valid_mega %>%
          filter(match == 1) %>%
          group_by(main_gene, main_cs_index) %>%
          filter(cos == max(cos)) %>%
          select(trait = main_gene, CSIndex = main_cs_index, match, cos, w_jac, n_jac))) 

sum(cov_mega_res$match, na.rm = TRUE)
sum(mega_res$match, na.rm = TRUE)

rep_comp_ratio <- rep_comp_ratio %>%
  bind_rows(
    tidy(prop.test(c(sum(cov_mega_res$match, na.rm = TRUE),
      sum(mega_res$match, na.rm = TRUE)),
      c(nrow(cov_mega_res), nrow(mega_res)))) %>%
      mutate(type = "Cov VS Mega"))

mean(cov_mega_res$cos, na.rm = TRUE)
mean(mega_res$cos, na.rm = TRUE)

rep_comp_cos <- rep_comp_cos %>%
  bind_rows(
    tidy(t.test(cov_mega_res$cos, mega_res$cos)) %>%
      mutate(type = "Cov VS Mega"))


df_comp <- rep_comp_ratio %>%
  select(type, estimate1, estimate2, p.value) %>%
  mutate(class = "z") %>%
  bind_rows(rep_comp_cos %>%
  select(type, estimate1, estimate2, p.value) %>%
      mutate(class = "t")) %>%
  mutate(type = factor(type, levels = c("Cov VS Indep", "Cov VS Meta", "Cov VS Mega")),
    class = factor(class, levels = c("z", "t"))) %>%
  arrange(type, class) %>%
  select(type, class, estimate1, estimate2, p.value)

write_tsv(df_comp, "./tables/s7.tsv")



