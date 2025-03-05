library(tidyverse)
library(glue)
library(broom)
library(ggpubr)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Documents/github/data/sushie_results/sim3"

# 2 pop general data
df_pip_pop2_tmp <- read_tsv(glue("{sim_data_path}/sushie_2pop_pip_all.tsv.gz"))

df_auprc <- tibble()
for (idx in seq(0.1, 0.9, 0.1))  {
  print(idx)
  df_tmp <- df_pip_pop2_tmp %>%
    select(-susie_ss_pip) %>%
    select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, causal, contains("pip")) %>%
    pivot_longer(cols = contains("pip")) %>%
    filter(!is.na(value)) %>%
    mutate(name = gsub("_pip", "", name)) %>%
    group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
    summarize(FP = sum(causal == 0 & value >= idx),
      TP = sum(causal == 1 & value >= idx),
      FN = sum(causal == 1 & value < idx),
      TN = sum(causal == 0 & value < idx),
      Precision = TP/(TP+FP),
      Recall = TP/(TP+FN)) %>%
    ungroup() %>%
    mutate(thresh = idx)
  df_auprc <- df_auprc %>%
    bind_rows(df_tmp)
}

pp_auprc <- df_auprc %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

# write_tsv(pp_auprc, "{sim_data_path}/auprc_data.tsv")

df_cali <- df_pip_pop2_tmp %>%
  select(-susie_ss_pip) %>%
  select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, causal, contains("pip")) %>%
  pivot_longer(cols = contains("pip")) %>%
  filter(!is.na(value)) %>%
  mutate(name = gsub("_pip", "", name)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  arrange(desc(value)) %>%
  mutate(bin = cut(value, breaks = seq(0, 1, length.out = 11),
    labels = FALSE, include.lowest = TRUE)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, bin) %>%
  summarize(m_val = mean(value),
    obs = mean(causal == 1),
    obs_se = sd(causal == 1)/sqrt(n())) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

# write_tsv(df_cali, "{sim_data_path}/calibration_data.tsv")


df_tmp1 <- df_pip_pop2_tmp %>%
  select(-susie_ss_pip) %>%
  select(sim, locus, N, L1, L2, L3, h2g, rho, causal, contains("pip")) %>%
  pivot_longer(cols = contains("pip")) %>%
  mutate(name = gsub("_pip", "", name)) %>%
  filter(!is.na(value)) %>%
  group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
  arrange(desc(value)) %>%
  mutate(
    False_Positives = cumsum(1 - value), 
    Total_Selected = row_number(),  
    Bayesian_FDR = False_Positives / Total_Selected
  )

optimal_threshold <- df_tmp1 %>%
  filter(Bayesian_FDR <= 0.05) %>%
  summarise(Optimal_PIP_Threshold = min(value, na.rm = TRUE))

df_pip_fdr <- df_tmp1 %>%
  left_join(optimal_threshold) %>%
  group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
  mutate(if_positive = value >= Optimal_PIP_Threshold & causal == 1,
    if_positive = ifelse(is.na(if_positive), FALSE, if_positive)) %>%
  summarize(value = as.numeric(sum(if_positive) > 0))

# write_tsv(df_pip_fdr, "{sim_data_path}/fdr_pip_data.tsv")

df_cs_pop2_tmp <- read_tsv(glue("{sim_data_path}/sushie_2pop_cs_all.tsv.gz"))

df_cs_fdr <- df_cs_pop2_tmp %>%
  filter(!method %in% "sushie_ss") %>%
  filter(!is.na(SNPIndex_1based)) %>%
  bind_rows(
    df_cs_pop2_tmp %>%
      filter(!method %in% "sushie_ss") %>%
      filter(is.na(SNPIndex_1based)) %>%
      distinct(SNPIndex_1based, method, sim, locus, N, L1, L2, L3, h2g, rho) %>%
      left_join(crossing(CSIndex = c(1, 2),
        method = c("sushie", "indep", "meta", "susie",
          "susiex", "mesusie",  "xmap",  "xmap_ind")),
        by = "method")
  ) %>%
  group_by(sim, locus, N, L1, L2, L3, h2g, rho, method) %>%
  summarize(Power = ifelse(is.na(SNPIndex_1based), 0, as.numeric(sum(causal == 1) != 0)))

# write_tsv(df_cs_fdr, "{sim_data_path}/fdr_cs_data.tsv")
