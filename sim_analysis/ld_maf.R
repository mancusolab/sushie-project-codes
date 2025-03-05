library(tidyverse)
library(broom)

# 2 pop general data
load("./data/df_2pop.RData")

# to replicate our analysis you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Downloads/sushie_sim_data_results"

# gene-level LD 
tmp_df_gene_ld <- read_tsv(glue("{sim_data_path}/sim_gene_ld.tsv.gz"))

df_gene_ld2 <- tmp_df_gene_ld %>%
  # second smalelst p values if p value is 0
  mutate(anova_p = ifelse(anova_p == 0, -log(8.683539e-213), -log(anova_p)),
    levene_p = -log(levene_p)) %>%
  select(locus = `0`, contains("stats"), contains("_p"), l2_diff, l2_diff_sq)

cor(df_gene_ld2[c(2:7)])
tidy(cor.test(df_gene_ld2$levene_p, df_gene_ld2$l2_diff))

tidy(cor.test(df_gene_ld2$anova_stats, df_gene_ld2$l2_diff))
tidy(cor.test(df_gene_ld2$anova_stats, df_gene_ld2$levene_p))

df_gene_ld <- tmp_df_gene_ld %>%
  # second smalelst p values if p value is 0
  mutate(anova_p = ifelse(anova_p == 0, -log(8.683539e-213), -log(anova_p)),
    levene_p = -log(levene_p)) %>%
  select(locus = `0`, contains("stats"), contains("_p"), l2_diff, l2_diff_sq) %>%
  pivot_longer(cols=-locus)

res_gene_ld <- tibble()
for (method_name in c("SuShiE", "SuShiE-Indep", "Meta-SuSiE",
  "SuSiE", "SuSiEx", "MESuSiE", "XMAP", "XMAP-IND")) {
  for (metric_name in unique(df_gene_ld$name)) {
    df_tmp <- df_2pop %>%
      filter(N %in% c("400:400")) %>%
      filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
      filter(h2g %in% "0.05:0.05" & rho == "0.8") %>%
      filter(name %in% method_name) %>%
      left_join(df_gene_ld %>%
          filter(name == metric_name) %>%
          rename(metric = name, 
            value2 = value), by = "locus")
    
    res_gene_ld <- res_gene_ld %>%
      bind_rows(
        tidy(lm(value ~ value2 +CSIndex ,
          filter(df_tmp, type == "PIP"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "PIP",
            name = method_name,
            metric = metric_name),
        tidy(lm(value ~ value2+CSIndex ,
          filter(df_tmp, type == "CS"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "CS",
            name = method_name,
            metric = metric_name),
        tidy(lm(value ~ value2+CSIndex ,
          filter(df_tmp, type == "Calibration"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "Calibration",
            name = method_name,
            metric = metric_name)
      ) 
  }
}

res_gene_ld %>%
  filter(metric %in% "levene_stats") %>%
  filter(name %in% "SuShiE")
0.0175/2
res_gene_ld %>%
  filter(metric %in% "l2_diff") %>%
  filter(name %in% "SuShiE")

res_gene_ld %>%
  filter(metric %in% "l2_diff_sq") %>%
  filter(name %in% "SuShiE")

res_gene_ld %>%
  filter(metric %in% "l2_diff")

res_gene_ld %>%
  filter(metric %in% "anova_stats") %>%
  filter(name %in% "SuShiE")
0.0871/2

df_causal_fst <- read_tsv(glue("{sim_data_path}/sim_causal_fst.tsv.gz"),
  col_names = FALSE) %>%
  select(locus = X3, sim = X2, fst = X1)

res_causal_fst <- tibble()
for (method_name in c("SuShiE", "SuShiE-Indep", "Meta-SuSiE",
  "SuSiE", "SuSiEx", "MESuSiE", "XMAP", "XMAP-IND")) {
  for (metric_name in unique(df_gene_ld$name)) {
    df_tmp <- df_2pop %>%
      filter(N %in% c("400:400")) %>%
      filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
      filter(h2g %in% "0.05:0.05" & rho == "0.8") %>%
      left_join(df_causal_fst, by = c( "sim", "locus"))
    
    fst_group <- df_tmp %>%
      distinct(sim, locus, fst) 
    
    low_fst_group <- fst_group %>%
      filter(fst >= 0.15) %>%
      pull(locus)
    
    high_fst_group <- fst_group %>%
      filter(fst <=0.05) %>%
      pull(locus)
    
    df_tmp1 <- df_2pop %>%
      filter(N %in% c("400:400")) %>%
      filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
      filter(h2g %in% "0.05:0.05" & rho == "0.8") %>%
      filter(name %in% method_name) %>%
      left_join(df_gene_ld %>%
          filter(name == metric_name) %>%
          rename(metric = name, 
            value2 = value), by = "locus") %>%
      filter(locus %in% low_fst_group)
    
    df_tmp2 <- df_2pop %>%
      filter(N %in% c("400:400")) %>%
      filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
      filter(h2g %in% "0.05:0.05" & rho == "0.8") %>%
      filter(name %in% method_name) %>%
      left_join(df_gene_ld %>%
          filter(name == metric_name) %>%
          rename(metric = name, 
            value2 = value), by = "locus") %>%
      filter(locus %in% high_fst_group)
    
    res_causal_fst <- res_causal_fst %>%
      bind_rows(
        tidy(lm(value ~ value2 +CSIndex,
          filter(df_tmp1, type == "PIP"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "PIP",
            name = method_name,
            metric = metric_name,
            fst_group = "low"),
        tidy(lm(value ~ value2+CSIndex,
          filter(df_tmp1, type == "CS"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "CS",
            name = method_name,
            metric = metric_name,
            fst_group = "low"),
        tidy(lm(value ~ value2+CSIndex,
          filter(df_tmp1, type == "Calibration"))) %>%
          filter(grepl("value", term)) %>%
          mutate(type = "Calibration",
            name = method_name,
            metric = metric_name,
            fst_group = "low"),
        tidy(lm(value ~ value2 +CSIndex,
          filter(df_tmp2, type == "PIP"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "PIP",
            name = method_name,
            metric = metric_name,
            fst_group = "high"),
        tidy(lm(value ~ value2+CSIndex,
          filter(df_tmp2, type == "CS"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "CS",
            name = method_name,
            metric = metric_name,
            fst_group = "high"),
        tidy(lm(value ~ value2+CSIndex,
          filter(df_tmp2, type == "Calibration"))) %>%
          filter(grepl("value2", term)) %>%
          mutate(type = "Calibration",
            name = method_name,
            metric = metric_name,
            fst_group = "high")
      )  
  }
}

res_causal_fst %>%
  filter(metric %in% "anova_stats") %>%
  filter(name %in% "SuShiE")

df_param <- df_2pop %>%
  distinct(sim, N, L1, L2, L3, h2g, rho)

tmp_df_fst <- read_tsv(glue("{sim_data_path}/sushie_real_fst.tsv.gz")) %>%
  filter(`#POP1` == "AFR" & POP2 == "EUR") %>%
  mutate(method = "Real Data from TOPMed MESA mRNA") %>%
  select(method, fst = `HUDSON_FST`) %>%
  filter(fst >=0)

set.seed(123)
df_fst <- tmp_df_fst %>%
  bind_rows(df_causal_fst %>%
      left_join(df_param) %>%
      sample_n(nrow(tmp_df_fst)) %>%
      mutate(method = "Sim Data from 1000G") %>%
      select(method, fst))

ggplot(df_fst, aes(x = fst)) +
  geom_histogram(binwidth = 0.01, fill = "blue", alpha = 0.6, color = "black") +
  labs(title = "Distribution of Fst Values",
    x = "Fst",
    y = "Frequency") +
  facet_wrap(~method) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = legend_fontsize),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=fontsize),
    axis.text=element_text(size = fontsize),
    legend.position = "bottom")

# ggsave("./manuscript_plots/additional/s10.png", width = 6, height = 4)

