# this script is to output numbers in the manuscript and 
# the first two supp tables

library(tidyverse)
library(broom)

# abstract and introduction
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.sushie_cs.tsv.gz")
length(unique(rnaseq_cov$trait))

proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.sushie_cs.tsv.gz")
length(unique(proteins_cov$trait))

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.sushie_cs.tsv.gz")
length(unique(genoa_cov$trait))

length(unique(rnaseq_cov$trait)) + length(unique(proteins_cov$trait)) + length(unique(genoa_cov$trait))

# validation
v5_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_normal.sushie_cs.tsv.gz")
length(unique(v5_cov$trait))

interval_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/interval_normal.sushie_cs.tsv.gz")
length(unique(interval_cov$trait))

geuvadis_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_normal.sushie_cs.tsv.gz")
length(unique(geuvadis_cov$trait))

length(unique(v5_cov$trait)) + length(unique(interval_cov$trait)) +
  length(unique(geuvadis_cov$trait))

# overlap
v5_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata2/v5_overlap_gene_list_noMHC.tsv", col_names = FALSE)

interval_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata2/interval_overlap_gene_list_noMHC.tsv", col_names = FALSE)

geuvadis_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata2/geuvadis_overlap_gene_list_noMHC.tsv", col_names = FALSE)


nrow(v5_meta %>% filter(X1 %in% v5_cov$trait))
nrow(interval_meta %>% filter(X2 %in% interval_cov$trait))
nrow(geuvadis_meta %>% filter(X1 %in% geuvadis_cov$trait))

nrow(v5_meta %>% filter(X1 %in% v5_cov$trait)) +
  nrow(interval_meta %>% filter(X2 %in% interval_cov$trait)) +
  nrow(geuvadis_meta %>% filter(X1 %in% geuvadis_cov$trait))


# more molQTLs
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.sushie_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_indep.sushie_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.meta_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = meta_pip_all)

rnaseq_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.mega_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_susiex_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_mesusie_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

# rnaseq_multisusie <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_multisusie_cs.tsv.gz") %>%
#   mutate(study = "mesa.mrna") %>%
#   select(study, trait, CSIndex, snp, pip)

rnaseq_xmap <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_xmap_cs.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.sushie_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

proteins_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_indep.sushie_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

proteins_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.meta_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = meta_pip_all)

proteins_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.mega_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

proteins_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_susiex_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

proteins_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_mesusie_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)
# 
# proteins_multisusie <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_multisusie_cs.tsv.gz") %>%
#   mutate(study = "mesa.protens") %>%
#   select(study, trait, CSIndex, snp, pip)

proteins_xmap <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_xmap_cs.tsv.gz") %>%
  mutate(study = "mesa.protens") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.sushie_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

genoa_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_indep.sushie_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

genoa_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.meta_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = meta_pip_all)

genoa_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.mega_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

genoa_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_susiex_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

genoa_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_mesusie_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

# genoa_multisusie <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_multisusie_cs.tsv.gz") %>%
#   mutate(study = "genoa.mrna") %>%
#   select(study, trait, CSIndex, snp, pip)

genoa_xmap <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_xmap_cs.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

cp0 <- nrow(distinct(filter(rnaseq_cov, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_cov, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_cov, !is.na(snp)), trait))

cp1 <- nrow(distinct(filter(rnaseq_indep, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_indep, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_indep, !is.na(snp)), trait))

cp2 <- nrow(distinct(filter(rnaseq_meta, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_meta, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_meta, !is.na(snp)), trait))

cp3 <- nrow(distinct(filter(rnaseq_mega, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_mega, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_mega, !is.na(snp)), trait))

cp4 <- nrow(distinct(filter(rnaseq_susiex, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_susiex, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_susiex, !is.na(snp)), trait))

cp5 <- nrow(distinct(filter(rnaseq_mesusie, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_mesusie, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_mesusie, !is.na(snp)), trait))

# cp6 <- nrow(distinct(filter(rnaseq_multisusie, !is.na(snp)), trait)) +
#   nrow(distinct(filter(genoa_multisusie, !is.na(snp)), trait)) +
#   nrow(distinct(filter(proteins_multisusie, !is.na(snp)), trait))

cp7 <- nrow(distinct(filter(rnaseq_xmap, !is.na(snp)), trait)) +
  nrow(distinct(filter(genoa_xmap, !is.na(snp)), trait)) +
  nrow(distinct(filter(proteins_xmap, !is.na(snp)), trait))


cp7

((cp0+cp1+cp2+cp3+cp4+cp5)/6)/(length(unique(rnaseq_cov$trait)) + length(unique(proteins_cov$trait)) +
    length(unique(genoa_cov$trait)))

cp7/(length(unique(rnaseq_cov$trait)) + length(unique(proteins_cov$trait)) +
    length(unique(genoa_cov$trait)))

cp0 

cp0 - (cp1 + cp2 + cp3 + cp4 +cp5)/5
cp0 / ( (cp1 + cp2 + cp3 + cp4 +cp5)/5)


n_cov <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

n_indep <- rnaseq_indep %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

n_meta <- rnaseq_meta %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

n_mega <- rnaseq_mega %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

n_susiex <- rnaseq_susiex %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

n_mesusie <- rnaseq_mesusie %>%
  filter(!is.na(snp)) %>%
  distinct(trait, CSIndex) %>%
  nrow()

# n_multisusie <- rnaseq_multisusie %>%
#   filter(!is.na(snp)) %>%
#   distinct(trait, CSIndex) %>%
#   nrow()

n_cov - (n_indep + n_meta + n_mega + n_susiex + n_mesusie)/5

n_cov /((n_indep + n_meta + n_mega + n_susiex + n_mesusie)/5)

sushie_comp <- rnaseq_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(proteins_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count),
    genoa_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count)

df_count_indep <- sushie_comp %>%
  full_join(rnaseq_indep %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(proteins_indep %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count),
        genoa_indep %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(indep = count),
    by = "trait")

df_count_mega <- sushie_comp %>%
  full_join(rnaseq_mega %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(proteins_mega %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count),
        genoa_mega %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(mega = count),
    by = "trait")

df_count_meta <- sushie_comp %>%
  full_join(rnaseq_meta %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(proteins_meta %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count),
        genoa_meta %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(meta = count),
    by = "trait")

rnaseq_all_traits <- distinct(rnaseq_cov, trait)
proteins_all_traits <- distinct(proteins_cov, trait)
genoa_all_traits <- distinct(genoa_cov, trait)

df_count_susiex <- sushie_comp %>%
  full_join(rnaseq_susiex %>%
      distinct(trait) %>%
      mutate(count = 1) %>%
      right_join(rnaseq_all_traits, by = "trait") %>%
      mutate(count = ifelse(is.na(count), 0, 1)) %>%
      bind_rows(
        proteins_susiex %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(proteins_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
        genoa_susiex %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(genoa_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
      ) %>%
      rename(susiex = count),
    by = "trait")

df_count_mesusie <- sushie_comp %>%
  full_join(rnaseq_mesusie %>%
      distinct(trait) %>%
      mutate(count = 1) %>%
      right_join(rnaseq_all_traits, by = "trait") %>%
      mutate(count = ifelse(is.na(count), 0, 1)) %>%
      bind_rows(
        proteins_mesusie %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(proteins_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
        genoa_mesusie %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(genoa_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
      ) %>%
      rename(mesusie = count),
    by = "trait")

# df_count_multisusie <- sushie_comp %>%
#   full_join(rnaseq_multisusie %>%
#       distinct(trait) %>%
#       mutate(count = 1) %>%
#       right_join(rnaseq_all_traits, by = "trait") %>%
#       mutate(count = ifelse(is.na(count), 0, 1)) %>%
#       bind_rows(
#         proteins_multisusie %>%
#           distinct(trait) %>%
#           mutate(count = 1) %>%
#           right_join(proteins_all_traits, by = "trait") %>%
#           mutate(count = ifelse(is.na(count), 0, 1)),
#         genoa_multisusie %>%
#           distinct(trait) %>%
#           mutate(count = 1) %>%
#           right_join(genoa_all_traits, by = "trait") %>%
#           mutate(count = ifelse(is.na(count), 0, 1)),
#       ) %>%
#       rename(multisusie = count),
#     by = "trait")

tidy(t.test(df_count_indep$sushie, df_count_indep$indep))
tidy(t.test(df_count_mega$sushie, df_count_mega$mega))
tidy(t.test(df_count_meta$sushie, df_count_meta$meta))
tidy(t.test(df_count_susiex$sushie, df_count_susiex$susiex))
tidy(t.test(df_count_mesusie$sushie, df_count_mesusie$mesusie))
# tidy(t.test(df_count_multisusie$sushie, df_count_multisusie$multisusie, alternative = "greater"))
cp0 / cp2



# basic summary data comparison
source("utils.R")
method_col <- tibble("method" = c("SuShiE-Cov", "SuShiE-Indep",
  "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE"))

comp_col <- tibble("type" = 1:5,
  comp = c("vs. SuShiE-Indep", "vs. Meta-SuSiE", "vs. SuSiE",
    "vs. SuSiEx", "vs. MESuSiE"))

n_total_rnaseq <- length(unique(rnaseq_cov$trait))
n_total_proteins <- length(unique(proteins_cov$trait))
n_total_genoa <- length(unique(genoa_cov$trait))
n_total <- n_total_rnaseq + n_total_proteins + n_total_genoa

# aggregate
agg_res <- round(basic_sum(
  bind_rows(rnaseq_cov, proteins_cov, genoa_cov),
  bind_rows(rnaseq_indep, proteins_indep, genoa_indep),
  bind_rows(rnaseq_meta, proteins_meta, genoa_meta),
  bind_rows(rnaseq_mega, proteins_mega, genoa_mega),
  bind_rows(rnaseq_susiex, proteins_susiex, genoa_susiex),
  bind_rows(rnaseq_mesusie, proteins_mesusie, genoa_mesusie),
  # bind_rows(rnaseq_multisusie, proteins_multisusie, genoa_multisusie),
  total = n_total
), 2) %>%
  mutate(type = row_number() - 1)

agg_comp <- basic_compare(
  bind_rows(rnaseq_cov, proteins_cov, genoa_cov),
  bind_rows(rnaseq_indep, proteins_indep, genoa_indep),
  bind_rows(rnaseq_meta, proteins_meta, genoa_meta),
  bind_rows(rnaseq_mega, proteins_mega, genoa_mega),
  bind_rows(rnaseq_susiex, proteins_susiex, genoa_susiex),
  bind_rows(rnaseq_mesusie, proteins_mesusie, genoa_mesusie)
  # bind_rows(rnaseq_multisusie, proteins_multisusie, genoa_multisusie)
  ) 

comp_stats1 <- agg_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(agg_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())


# rnaseq
rnaseq_res <- round(basic_sum(rnaseq_cov, rnaseq_indep,
  rnaseq_meta, rnaseq_mega, rnaseq_susiex, rnaseq_mesusie,
   total=n_total_rnaseq), 2) %>%
  mutate(type = row_number() - 1)

rnaseq_comp <- basic_compare(rnaseq_cov, rnaseq_indep, rnaseq_meta, rnaseq_mega,
  rnaseq_susiex, rnaseq_mesusie) 

comp_stats2 <- rnaseq_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(rnaseq_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())

# proteins
proteins_res <- round(basic_sum(proteins_cov, proteins_indep,
  proteins_meta, proteins_mega, proteins_susiex,
  proteins_mesusie, total = n_total_proteins), 2) %>%
  mutate(type = row_number() - 1)

proteins_comp <- basic_compare(proteins_cov,
  proteins_indep, proteins_meta, proteins_mega,
  proteins_susiex, proteins_mesusie) 

comp_stats3 <- proteins_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(proteins_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())

# genoa
genoa_res <- round(basic_sum(genoa_cov, genoa_indep,
  genoa_meta, genoa_mega, genoa_susiex, genoa_mesusie,
  total = n_total_genoa), 2) %>%
  mutate(type = row_number() - 1)

genoa_comp <- basic_compare(genoa_cov, genoa_indep,
  genoa_meta, genoa_mega, genoa_susiex, genoa_mesusie) 

comp_stats4 <- genoa_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(genoa_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())

comp_stats <- comp_stats1 %>%
  left_join(comp_stats2, by = "name") %>%
  left_join(comp_stats3, by = "name") %>%
  left_join(comp_stats4, by = "name")

comp_stats %>%
  filter(grepl("estimate_cs_size", name)) %>%
  summarize(value = mean(value.x))

comp_stats %>%
  filter(grepl("estimate_avg_pip", name)) %>%
  summarize(value = mean(value.x))

comp_stats %>%
  filter(grepl("estimate_g09", name)) %>%
  summarize(value = mean(value.x))

# Define the desired order of row names
desired_order <- c(
  "cs_num",
  "n_cs_size",
  "estimate_cs_size",
  "p.value_cs_size",
  "estimate_avg_pip",
  "p.value_avg_pip",
  "estimate_g09",
  "p.value_g09"
)

df1 <- comp_stats[1:5,] %>%
  bind_rows(comp_stats %>%
      filter(row_number() > 5) %>%
      filter(!grepl("n_avg_pip", name)) %>%
      filter(!grepl("n_g09", name)) %>%
      separate(name, into = c("name", "number"), sep = "_(?=[^_]+$)", convert = TRUE) %>%
      filter(!name %in% c("cs_size", "avg_pip", "g09", "total")) %>%
      group_by(number) %>%
      mutate(name = factor(name, levels = desired_order)) %>%
      arrange(number, name) %>%
      ungroup() %>%
      select(-number),
  )

write_tsv(df1, "./tables/s4.tsv")

# % of 1-3 molQTLs 
rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  bind_rows(proteins_cov %>%
      filter(!is.na(snp)) %>%
      group_by(trait) %>%
      summarize(n = length(unique(CSIndex)))) %>%
  bind_rows(genoa_cov %>%
      filter(!is.na(snp)) %>%
      group_by(trait) %>%
      summarize(n = length(unique(CSIndex)))) %>%
  mutate(total = n()) %>%
  group_by(n) %>%
  summarize(n_total = n(),
    total = first(total)) %>%
  mutate(c_total = cumsum(n_total),
    prop = c_total/total) 


# study specific
nrow(distinct(filter(rnaseq_cov, !is.na(snp)), trait))
nrow(distinct(filter(proteins_cov, !is.na(snp)), trait))
nrow(distinct(filter(genoa_cov, !is.na(snp)), trait))


rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(total = n()) %>%
  group_by(n) %>%
  summarize(n_total = n(),
    total = first(total)) %>%
  mutate(c_total = cumsum(n_total),
    prop = c_total/total) 

proteins_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(total = n()) %>%
  group_by(n) %>%
  summarize(n_total = n(),
    total = first(total)) %>%
  mutate(c_total = cumsum(n_total),
    prop = c_total/total) 

genoa_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(total = n()) %>%
  group_by(n) %>%
  summarize(n_total = n(),
    total = first(total)) %>%
  mutate(c_total = cumsum(n_total),
    prop = c_total/total) 

# validation 
v5_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_normal.sushie_cs.tsv.gz")  %>%
  mutate(study = "mesa.v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

v5_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_indep.sushie_cs.tsv.gz")  %>%
  mutate(study = "mesa.v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

v5_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_normal.meta_cs.tsv.gz")  %>%
  mutate(study = "mesa.v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

v5_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_normal.mega_cs.tsv.gz")  %>%
  mutate(study = "mesa.v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

v5_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_susiex_cs.tsv.gz") %>%
  mutate(study = "v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

v5_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/v5_mesusie_cs.tsv.gz") %>%
  mutate(study = "v5") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)


geuvadis_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_normal.sushie_cs.tsv.gz")  %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

geuvadis_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_indep.sushie_cs.tsv.gz") %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

geuvadis_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_normal.meta_cs.tsv.gz") %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

geuvadis_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_normal.mega_cs.tsv.gz") %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

geuvadis_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_susiex_cs.tsv.gz") %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

geuvadis_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/geuvadis_mesusie_cs.tsv.gz") %>%
  mutate(study = "geuvadis.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)


cp0 <- nrow(distinct(filter(v5_cov, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_cov, !is.na(snp)), trait))

cp1 <- nrow(distinct(filter(v5_indep, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_indep, !is.na(snp)), trait)) 

cp2 <- nrow(distinct(filter(v5_meta, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_meta, !is.na(snp)), trait))

cp3 <- nrow(distinct(filter(v5_mega, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_mega, !is.na(snp)), trait))

cp4 <- nrow(distinct(v5_susiex, trait)) +
  nrow(distinct(geuvadis_susiex, trait))

cp5 <- nrow(distinct(v5_mesusie, trait)) +
  nrow(distinct(geuvadis_mesusie, trait))

cp0 - (cp1 + cp2 + cp3 + cp4 + cp5)/5
cp0 /( (cp1 + cp2 + cp3 + cp4 +cp5 )/5)

sushie_comp <- v5_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(
    geuvadis_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count)

df_count_indep <- sushie_comp %>%
  full_join(v5_indep %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(
        geuvadis_indep %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(indep = count),
    by = "trait")

df_count_mega <- sushie_comp %>%
  full_join(v5_mega %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(
        geuvadis_mega %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(mega = count),
    by = "trait")

df_count_meta <- sushie_comp %>%
  full_join(v5_meta %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(
        geuvadis_meta %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(meta = count),
    by = "trait")

v5_all_traits <- distinct(v5_cov, trait)
geuvadis_all_traits <- distinct(geuvadis_cov, trait)

df_count_susiex <- sushie_comp %>%
  full_join(v5_susiex %>%
      distinct(trait) %>%
      mutate(count = 1) %>%
      right_join(v5_all_traits, by = "trait") %>%
      mutate(count = ifelse(is.na(count), 0, 1)) %>%
      bind_rows(
        geuvadis_susiex %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(geuvadis_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
      ) %>%
      rename(susiex = count),
    by = "trait")

df_count_mesusie <- sushie_comp %>%
  full_join(v5_mesusie %>%
      distinct(trait) %>%
      mutate(count = 1) %>%
      right_join(v5_all_traits, by = "trait") %>%
      mutate(count = ifelse(is.na(count), 0, 1)) %>%
      bind_rows(
        geuvadis_mesusie %>%
          distinct(trait) %>%
          mutate(count = 1) %>%
          right_join(geuvadis_all_traits, by = "trait") %>%
          mutate(count = ifelse(is.na(count), 0, 1)),
      ) %>%
      rename(mesusie = count),
    by = "trait")

tidy(t.test(df_count_indep$sushie, df_count_indep$indep))
tidy(t.test(df_count_meta$sushie, df_count_meta$meta))
tidy(t.test(df_count_mega$sushie, df_count_mega$mega))
tidy(t.test(df_count_susiex$sushie, df_count_susiex$susiex))
tidy(t.test(df_count_mesusie$sushie, df_count_mesusie$mesusie))

method_col <- tibble("method" = c("SuShiE-Cov", "SuShiE-Indep",
  "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE"))

comp_col <- tibble("type" = 1:5,
  comp = c("vs. SuShiE-Indep", "vs. Meta-SuSiE", "vs. SuSiE",
    "vs. SuSiEx", "vs. MESuSiE"))

n_total_v5 <- length(unique(v5_cov$trait))
n_total_geuvadis <- length(unique(geuvadis_cov$trait))
n_total <- n_total_v5+ n_total_geuvadis

# aggregate
agg_res <- round(basic_sum(
  bind_rows(v5_cov,  geuvadis_cov),
  bind_rows(v5_indep,  geuvadis_indep),
  bind_rows(v5_meta, geuvadis_meta),
  bind_rows(v5_mega, geuvadis_mega),
  bind_rows(v5_susiex, geuvadis_susiex),
  bind_rows(v5_mesusie, geuvadis_mesusie),
  total = n_total
), 2) %>%
  mutate(type = row_number() - 1)

agg_comp <- basic_compare(
  bind_rows(v5_cov,  geuvadis_cov),
  bind_rows(v5_indep, geuvadis_indep),
  bind_rows(v5_meta, geuvadis_meta),
  bind_rows(v5_mega, geuvadis_mega),
  bind_rows(v5_susiex, geuvadis_susiex),
  bind_rows(v5_mesusie, geuvadis_mesusie)
) 


comp_stats1 <- agg_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(agg_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())


# v5
v5_res <- round(basic_sum(v5_cov, v5_indep,
  v5_meta, v5_mega, v5_susiex, v5_mesusie,
  total=n_total_v5), 2) %>%
  mutate(type = row_number() - 1)

v5_comp <- basic_compare(v5_cov, v5_indep, v5_meta, v5_mega,
  v5_susiex, v5_mesusie) 

comp_stats2 <- v5_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(v5_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())

# geuvadis
geuvadis_res <- round(basic_sum(geuvadis_cov, geuvadis_indep,
  geuvadis_meta, geuvadis_mega, geuvadis_susiex, geuvadis_mesusie,
  total = n_total_geuvadis), 2) %>%
  mutate(type = row_number() - 1)

geuvadis_comp <- basic_compare(geuvadis_cov, geuvadis_indep,
  geuvadis_meta, geuvadis_mega, geuvadis_susiex, geuvadis_mesusie) 

comp_stats3 <- geuvadis_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(geuvadis_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"),
        contains("4"), contains("5"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3"),
    contains("_4"), contains("_5")) %>%
  pivot_longer(cols=everything())

comp_stats <- comp_stats1 %>%
  left_join(comp_stats2, by = "name") %>%
  left_join(comp_stats3, by = "name")

comp_stats %>%
  filter(grepl("estimate_cs_size", name)) %>%
  summarize(value = mean(value.x))

comp_stats %>%
  filter(grepl("estimate_avg_pip", name)) %>%
  summarize(value = mean(value.x))

comp_stats %>%
  filter(grepl("estimate_g09", name)) %>%
  summarize(value = mean(value.x))

# Define the desired order of row names
desired_order <- c(
  "cs_num",
  "n_cs_size",
  "estimate_cs_size",
  "p.value_cs_size",
  "estimate_avg_pip",
  "p.value_avg_pip",
  "estimate_g09",
  "p.value_g09"
)

df1 <- comp_stats[1:5,] %>%
  bind_rows(comp_stats %>%
      filter(row_number() > 5) %>%
      filter(!grepl("n_avg_pip", name)) %>%
      filter(!grepl("n_g09", name)) %>%
      separate(name, into = c("name", "number"), sep = "_(?=[^_]+$)", convert = TRUE) %>%
      filter(!name %in% c("cs_size", "avg_pip", "g09", "total")) %>%
      group_by(number) %>%
      mutate(name = factor(name, levels = desired_order)) %>%
      arrange(number, name) %>%
      ungroup() %>%
      select(-number),
  )

# write_tsv(df1 %>% select(-name), "./tables/s8.tsv", col_names = FALSE)

interval_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/interval_normal.sushie_cs.tsv.gz") %>%
  mutate(study = "interval") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

n_total_interval <- length(unique(interval_cov$trait))

basic_sum(bind_rows(v5_cov, geuvadis_cov, interval_cov))
basic_sum(interval_cov, total = n_total_interval)

loss_ct <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_loss.tsv.gz")

loss_ct  %>%
  pivot_longer(cols = -1) %>%
  group_by(name) %>%
  summarize(mval = mean(value))

2566-2537
2537-2036
