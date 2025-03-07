library(tidyverse)
library(glue)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Downloads/sushie_sim_data_results/"

# 2 pop general data
sushie_pip_pop2 <- read_tsv(glue("{sim_data_path}/sushie_2pop_pip.tsv.gz"))

susiex_in_pip_pop2 <- read_tsv(glue("{sim_data_path}/susiex_in_2pop_pip.tsv.gz"))

mesusie_in_pip_pop2 <- read_tsv(glue("{sim_data_path}/mesusie_in_2pop_pip.tsv.gz"))

xmap_in_pip_pop2 <- read_tsv(glue("{sim_data_path}/xmap_in_2pop_pip.tsv.gz"))

xmap_ind_pip_pop2 <- read_tsv(glue("{sim_data_path}/xmap_ind_2pop_pip.tsv.gz"))

sushie_cs_pop2 <- read_tsv(glue("{sim_data_path}/sushie_2pop_cs.tsv.gz"))

susiex_in_cs_pop2 <- read_tsv(glue("{sim_data_path}/susiex_in_2pop_cs.tsv.gz"))

mesusie_in_cs_pop2 <- read_tsv(glue("{sim_data_path}/mesusie_in_2pop_cs.tsv.gz"))

xmap_in_cs_pop2 <- read_tsv(glue("{sim_data_path}/xmap_in_2pop_cs.tsv.gz"))

xmap_ind_cs_pop2 <- read_tsv(glue("{sim_data_path}/xmap_ind_2pop_cs.tsv.gz"))

ref_param <- read_tsv(glue("{sim_data_path}/sushie_2pop_pip.tsv.gz")) %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho)

dd_pip_pop2 <- sushie_pip_pop2 %>%
  select(sushie_pip, indep_pip, meta_pip, susie_pip,
    CSIndex, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = contains("pip")) %>%
  mutate(name = gsub("_pip", "", name)) %>%
  bind_rows(
    susiex_in_pip_pop2 %>%
      select(pip, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "susiex") %>%
      rename(value = pip),
    mesusie_in_pip_pop2 %>%
      select(mesusie, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "mesusie") %>%
      rename(value = mesusie),
    xmap_in_pip_pop2 %>%
      select(xmap, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap") %>%
      rename(value = xmap),
    xmap_ind_pip_pop2 %>%
      select(xmap, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap_ind") %>%
      rename(value = xmap)
  ) %>%
  filter(!is.na(value))

dd_cali_pop2 <- sushie_pip_pop2 %>%
  select(sushie_cali, indep_cali, meta_cali, susie_cali,
    sim, locus, N, L1, L2, L3, h2g, rho, CSIndex) %>%
  pivot_longer(cols = contains("cali")) %>%
  mutate(name = gsub("_cali", "", name)) %>%
  bind_rows(
    susiex_in_pip_pop2 %>%
      select(cali, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "susiex") %>%
      rename(value = cali),
    mesusie_in_pip_pop2 %>%
      select(cali, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "mesusie") %>%
      rename(value = cali),
    xmap_in_pip_pop2 %>%
      select(cali, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap") %>%
      rename(value = cali),
    xmap_ind_pip_pop2 %>%
      select(cali, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap_ind") %>%
      rename(value = cali)
  ) %>%
  filter(!is.na(value))

dd_cs_pop2 <- sushie_cs_pop2 %>%
  select(sushie:susie, sim, locus, N, L1, L2, L3, h2g, rho, CSIndex) %>%
  pivot_longer(cols = sushie:susie) %>%
  bind_rows(
    susiex_in_cs_pop2 %>%
      select(susiex, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "susiex") %>%
      rename(value = susiex),
    mesusie_in_cs_pop2 %>%
      select(mesusie, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "mesusie") %>%
      rename(value = mesusie),
    xmap_in_cs_pop2 %>%
      select(xmap, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap") %>%
      rename(value = xmap),
    xmap_ind_cs_pop2 %>%
      select(xmap, sim, locus, CSIndex) %>%
      left_join(ref_param) %>%
      mutate(name = "xmap_ind") %>%
      rename(value = xmap)
  ) %>%
  filter(!is.na(value))

df_2pop <- dd_pip_pop2 %>%
  mutate(type = "PIP") %>%
  bind_rows(dd_cs_pop2 %>%
      mutate(type = "CS"),
    dd_cali_pop2 %>%
      mutate(type = "Calibration")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

save(df_2pop, file = "./data/df_2pop.RData")

# 3 pop

pop3_pip <- read_tsv(glue("{sim_data_path}/sim_3pop_pip.tsv.gz"))
pop3_cs <- read_tsv(glue("{sim_data_path}/sim_3pop_cs.tsv.gz"))

dd_pip_pop3 <- pop3_pip %>%
  pivot_longer(cols = contains("pip")) %>%
  mutate(name = gsub("_pip", "", name)) %>%
  filter(!is.na(value)) %>%
  select(sim, locus, N, name, value) %>%
  mutate(type = "PIP")

dd_cs_pop3 <- pop3_cs %>%
  pivot_longer(cols = sushie1:sushie3) %>%
  filter(!is.na(value)) %>%
  select(sim, locus, N, name, value) %>%
  mutate(type = "CS")

dd_cali_pop3 <- pop3_pip %>%
  pivot_longer(cols = contains("cali")) %>%
  mutate(name = gsub("_cali", "", name)) %>%
  filter(!is.na(value)) %>%
  select(sim, locus, N, name, value) %>%
  mutate(type = "Calibration")

df_3pop <- bind_rows(dd_pip_pop3, dd_cs_pop3, dd_cali_pop3) %>%
  mutate(name = factor(name,
    labels = c("One-Ancestry", "Two-Ancestry", "Three-Ancestry"),
    levels = c("sushie1", "sushie2", "sushie3")))

save(df_3pop, file = "./data/df_3pop.RData")

# no shared
sushie_noshared <- read_tsv(glue("{sim_data_path}/sushie_noshared_fdr.tsv.gz"))

mesusie_noshared <- read_tsv(glue("{sim_data_path}/mesusie_in_noshared_fdr.tsv.gz"))

susiex_noshared <- read_tsv(glue("{sim_data_path}/susiex.in_noshared_fdr.tsv.gz"))

xmap_ind_noshared <- read_tsv(glue("{sim_data_path}/xmap_ind_noshared_fdr.tsv.gz"))

xmap_noshared <- read_tsv(glue("{sim_data_path}/xmap_in_noshared_fdr.tsv.gz"))

noshared_param <- sushie_noshared %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho)

df_noshared <- sushie_noshared %>%
  filter(method %in% c("sushie", "indep", "meta", "susie")) %>%
  mutate(value = as.numeric(fdr_cs - as_in > 0)) %>%
  select(sim, locus, N, L1, L2, L3, h2g, rho, name = method, value) %>%
  bind_rows(
    mesusie_noshared %>%
      mutate(value = as.numeric(fdr_cs - as_in > 0)) %>%
      select(sim, locus, name = method, value) %>%
      inner_join(noshared_param, by = c("sim", "locus")),
    susiex_noshared %>%
      mutate(value = as.numeric(fdr_cs - as_in > 0)) %>%
      select(sim, locus, name = method, value) %>%
      inner_join(noshared_param, by = c("sim", "locus")),
    xmap_noshared %>%
      mutate(value = as.numeric(fdr_cs - as_in > 0)) %>%
      select(sim, locus, name = method, value) %>%
      inner_join(noshared_param, by = c("sim", "locus")),
    xmap_ind_noshared %>%
      mutate(value = as.numeric(fdr_cs - as_in > 0),
        name = "xmap_ind") %>%
      select(sim, locus, name, value) %>%
      inner_join(noshared_param, by = c("sim", "locus"))
  ) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

save(df_noshared, file = "./data/df_noshared.RData")

# r2 and TWAS
tmp_r2 <- read_tsv(glue("{sim_data_path}/sim_pred_r2.tsv.gz"))
tmp_twas <- read_tsv(glue("{sim_data_path}/sim_pred_twas.tsv.gz"))

df_r2 <- tmp_r2 %>%
  filter(method %in% c("sushie", "indep", "meta",
    "susie",  "mesusie.in", "xmap.in", "xmap.ind",
    "enet", "lasso", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  filter(!is.na(value)) %>%
  mutate(ancestry = ifelse(name %in% "ancestry1_weight1", 1, 2)) %>%
  select(sim, locus, N, ngwas, h2g, h2ge, ancestry, name = method, value) %>%
  mutate(type = "R2")

df_twas <- tmp_twas %>%
  select(sim, locus, N, ngwas, h2g, h2ge, ancestry, sushie, indep, meta, susie,
    mesusie.in, xmap.in, xmap.ind, enet, lasso, ridge) %>%
  pivot_longer(cols=c(sushie:ridge)) %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, ngwas, h2g, h2ge, name, ancestry) %>%
  mutate(pval = 2 * pnorm(abs(value), lower.tail = FALSE),
    adjpval = p.adjust(pval),
    n = n(),
    value = as.numeric(adjpval < 0.05)) %>% 
  select(sim, locus, N, ngwas, h2g, h2ge, ancestry, name, value) %>%
  mutate(type = "TWAS")

df_r2twas <- bind_rows(df_r2, df_twas) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "mesusie.in", "xmap.in", "xmap.ind",
      "enet", "lasso", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "MESuSiE", "XMAP", "XMAP-IND",
      "Elastic Net", "LASSO", "gBLUP")))

save(df_r2twas, file = "./data/df_r2twas.RData")

# SuShiE vs SuShiE-SS

df_sushie_comp <- sushie_pip_pop2 %>%
  select(sushie = sushie_pip, sushie_ss = susie_ss_pip,
    CSIndex, sim, locus, N, L1, L2, L3, h2g, rho)

save(df_sushie_comp, file = "./data/df_sushie_comp.RData")

