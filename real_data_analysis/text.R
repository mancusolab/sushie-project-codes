# this script is to output numbers in the manuscript and 
# the first two supp tables

library(tidyverse)
library(broom)

# abstract and introduction
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz")
length(unique(rnaseq_cov$trait))

proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz")
length(unique(proteins_cov$trait))

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz")
length(unique(genoa_cov$trait))

length(unique(rnaseq_cov$trait)) + length(unique(proteins_cov$trait)) + length(unique(genoa_cov$trait))

# validation
v5_cov <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.sushie_cs.tsv.gz")
length(unique(v5_cov$trait))

interval_cov <- read_tsv("~/Documents/github/data/sushie_results/real/interval_normal.sushie_cs.tsv.gz")
length(unique(interval_cov$trait))

geuvadis_cov <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.sushie_cs.tsv.gz")
length(unique(geuvadis_cov$trait))

length(unique(v5_cov$trait)) + length(unique(interval_cov$trait)) +
  length(unique(geuvadis_cov$trait))  

# overlap
v5_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata/v5_overlap_gene_list_noMHC.tsv", col_names = FALSE)

interval_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata/interval_overlap_gene_list_noMHC.tsv", col_names = FALSE)

geuvadis_meta <- read_tsv("~/Documents/github/data/sushie_results/metadata/geuvadis_overlap_gene_list_noMHC.tsv", col_names = FALSE)


nrow(v5_meta %>% filter(X1 %in% v5_cov$trait))
nrow(interval_meta %>% filter(X2 %in% interval_cov$trait))
nrow(geuvadis_meta %>% filter(X1 %in% geuvadis_cov$trait))

nrow(v5_meta %>% filter(X1 %in% v5_cov$trait)) +
nrow(interval_meta %>% filter(X2 %in% interval_cov$trait)) +
nrow(geuvadis_meta %>% filter(X1 %in% geuvadis_cov$trait))


# more molQTLs
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz")

rnaseq_indep <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_indep.sushie_cs.tsv.gz")

rnaseq_meta <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.meta_cs.tsv.gz")

rnaseq_mega <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.mega_cs.tsv.gz")


proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz")

proteins_indep <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_indep.sushie_cs.tsv.gz")

proteins_meta <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.meta_cs.tsv.gz")

proteins_mega <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.mega_cs.tsv.gz")


genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz")

genoa_indep <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_indep.sushie_cs.tsv.gz")

genoa_meta <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.meta_cs.tsv.gz")

genoa_mega <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.mega_cs.tsv.gz")


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

cp0 - (cp1 + cp2 + cp3)/3
(cp0 - (cp1 + cp2 + cp3)/3)/cp0

df_count_indep <- rnaseq_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(proteins_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count),
    genoa_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
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

tidy(t.test(df_count_indep$sushie, df_count_indep$indep, alternative = "greater"))

df_count_mega <- rnaseq_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(proteins_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count),
    genoa_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
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

tidy(t.test(df_count_mega$sushie, df_count_mega$mega, alternative = "greater"))

(cp0 - cp2)/cp0

df_count <- rnaseq_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(proteins_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count),
    genoa_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
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

tidy(t.test(df_count$sushie, df_count$meta, alternative = "greater"))


# basic summary data comparison
source("utils.R")
method_col <- tibble("method" = c("SuShiE-Cov", "SuShiE-Indep",
  "Meta-SuSiE", "SuSiE"))

comp_col <- tibble("type" = 1:3,
  comp = c("vs. SuShiE-Indep", "vs. Meta-SuSiE", "vs. SuSiE"))

base <- bind_rows(rnaseq_cov %>%
    mutate(study = "mesa.mrna"),
  proteins_cov %>%
    mutate(study = "mesa.protens"),
  genoa_cov %>%
    mutate(study = "genoa.mrna")) 

dd <- bind_rows(rnaseq_mega %>%
    mutate(study = "mesa.mrna"),
  proteins_mega %>%
    mutate(study = "mesa.protens"),
  genoa_mega %>%
    mutate(study = "genoa.mrna")) 


# aggregate
agg_res <- round(basic_sum(
  bind_rows(rnaseq_cov, proteins_cov, genoa_cov),
  bind_rows(rnaseq_indep, proteins_indep, genoa_indep),
  bind_rows(rnaseq_meta, proteins_meta, genoa_meta),
  bind_rows(rnaseq_mega, proteins_mega, genoa_mega)
), 2) %>%
  mutate(type = row_number() - 1)

agg_comp <- basic_compare(
  bind_rows(rnaseq_cov, proteins_cov, genoa_cov),
  bind_rows(rnaseq_indep, proteins_indep, genoa_indep),
  bind_rows(rnaseq_meta, proteins_meta, genoa_meta),
  bind_rows(rnaseq_mega, proteins_mega, genoa_mega)) 

comp_stats1 <- agg_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(agg_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


# rnaseq
rnaseq_res <- round(basic_sum(rnaseq_cov, rnaseq_indep,
  rnaseq_meta, rnaseq_mega), 2) %>%
  mutate(type = row_number() - 1)

rnaseq_comp <- basic_compare(rnaseq_cov, rnaseq_indep, rnaseq_meta, rnaseq_mega) 

comp_stats2 <- rnaseq_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(rnaseq_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


# proteins
proteins_res <- round(basic_sum(proteins_cov, proteins_indep,
  proteins_meta, proteins_mega), 2) %>%
  mutate(type = row_number() - 1)

proteins_comp <- basic_compare(proteins_cov,
  proteins_indep, proteins_meta, proteins_mega) 

comp_stats3 <- proteins_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(proteins_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


# genoa
genoa_res <- round(basic_sum(genoa_cov, genoa_indep,
  genoa_meta, genoa_mega), 2) %>%
  mutate(type = row_number() - 1)

genoa_comp <- basic_compare(genoa_cov, genoa_indep,
  genoa_meta, genoa_mega) 

comp_stats4 <- genoa_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(genoa_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())

comp_stats <- comp_stats1 %>%
  left_join(comp_stats2, by = "name") %>%
  left_join(comp_stats3, by = "name") %>%
  left_join(comp_stats4, by = "name")

write_tsv(comp_stats, "./tables/s2.tsv")

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
nrow(distinct(filter(genoa_cov, !is.na(snp)), trait))
nrow(distinct(filter(proteins_cov, !is.na(snp)), trait))

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
v5_cov <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.sushie_cs.tsv.gz")

v5_indep <- read_tsv("~/Documents/github/data/sushie_results/real/v5_indep.sushie_cs.tsv.gz")

v5_meta <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.meta_cs.tsv.gz")

v5_mega <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.mega_cs.tsv.gz")


geuvadis_cov <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.sushie_cs.tsv.gz")

geuvadis_indep <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_indep.sushie_cs.tsv.gz")

geuvadis_meta <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.meta_cs.tsv.gz")

geuvadis_mega <- read_tsv("~/Documents/github/data/sushie_results/real/geuvadis_normal.mega_cs.tsv.gz")

cp0 <- nrow(distinct(filter(v5_cov, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_cov, !is.na(snp)), trait))

cp1 <- nrow(distinct(filter(v5_indep, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_indep, !is.na(snp)), trait)) 

cp2 <- nrow(distinct(filter(v5_meta, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_meta, !is.na(snp)), trait))

cp3 <- nrow(distinct(filter(v5_mega, !is.na(snp)), trait)) +
  nrow(distinct(filter(geuvadis_mega, !is.na(snp)), trait))

cp0 - (cp1 + cp2 + cp3)/3
(cp0 - (cp1 + cp2 + cp3)/3)/cp0


df_count_indep <- v5_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(geuvadis_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
  full_join(v5_indep %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(geuvadis_indep %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(indep = count),
    by = "trait")

tidy(t.test(df_count_indep$sushie, df_count_indep$indep, alternative = "greater"))

df_count_meta <- v5_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(geuvadis_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
  full_join(v5_meta %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(geuvadis_meta %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(meta = count),
    by = "trait")

tidy(t.test(df_count_meta$sushie, df_count_meta$meta, alternative = "greater"))

df_count_mega <- v5_cov %>%
  mutate(count = ifelse(is.na(snp),0, 1)) %>%
  distinct(trait, count) %>%
  bind_rows(geuvadis_cov %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count)) %>%
  rename(sushie = count) %>%
  full_join(v5_mega %>%
      mutate(count = ifelse(is.na(snp),0, 1)) %>%
      distinct(trait, count) %>%
      bind_rows(geuvadis_mega %>%
          mutate(count = ifelse(is.na(snp),0, 1)) %>%
          distinct(trait, count)) %>%
      rename(mega = count),
    by = "trait")

tidy(t.test(df_count_mega$sushie, df_count_mega$mega, alternative = "greater"))


method_col <- tibble("method" = c("SuShiE-Cov", "SuShiE-Indep",
  "Meta-SuSiE", "SuSiE"))

comp_col <- tibble("type" = 1:3,
  comp = c("vs. SuShiE-Indep", "vs. Meta-SuSiE", "vs. SuSiE"))

# aggregate
agg_res <- round(basic_sum(
  bind_rows(v5_cov, geuvadis_cov),
  bind_rows(v5_indep, geuvadis_indep),
  bind_rows(v5_meta, geuvadis_meta),
  bind_rows(v5_mega, geuvadis_mega)
), 2) %>%
  mutate(type = row_number() - 1)

agg_comp <- basic_compare(
  bind_rows(v5_cov, geuvadis_cov),
  bind_rows(v5_indep, geuvadis_indep),
  bind_rows(v5_meta, geuvadis_meta),
  bind_rows(v5_mega, geuvadis_mega)
)

comp_stats1 <- agg_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(agg_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


# v5
v5_res <- round(basic_sum(v5_cov, v5_indep,
  v5_meta, v5_mega), 2) %>%
  mutate(type = row_number() - 1)

v5_comp <- basic_compare(v5_cov, v5_indep, v5_meta, v5_mega) 

comp_stats2 <- v5_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(v5_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


# proteins
geuvadis_res <- round(basic_sum(geuvadis_cov, geuvadis_indep,
  geuvadis_meta, geuvadis_mega), 2) %>%
  mutate(type = row_number() - 1)

geuvadis_comp <- basic_compare(geuvadis_cov,
  geuvadis_indep, geuvadis_meta, geuvadis_mega) 

comp_stats3 <- geuvadis_res %>%
  pivot_wider(names_from = c(type),
    values_from = c(total, cs_num, cs_size, avg_pip, g09)) %>%
  bind_cols(geuvadis_comp %>%
      pivot_wider(names_from = c(metric, type),
        values_from = c(estimate, p.value, n)) %>%
      select(contains("1"), contains("2"), contains("3"))) %>%
  select(contains("_0"), contains("_1"), contains("_2"), contains("_3")) %>%
  pivot_longer(cols=everything())


comp_stats <- comp_stats1 %>%
  left_join(comp_stats2, by = "name") %>%
  left_join(comp_stats3, by = "name")

write_tsv(comp_stats, "./tables/s6.tsv")

interval_cov <- read_tsv("~/Documents/github/data/sushie_results/real/interval_normal.sushie_cs.tsv.gz")

basic_sum(bind_rows(v5_cov, geuvadis_cov, interval_cov))
basic_sum(interval_cov)

loss_ct <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_loss.tsv.gz")

loss_ct  %>%
  pivot_longer(cols = -1) %>%
  group_by(name) %>%
  summarize(mval = mean(value))

2566-2537
2537-2036
