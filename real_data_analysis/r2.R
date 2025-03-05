library(tidyverse)
library(ggpubr)
library(broom)
source("./utils.R")

# change the data folder to the zenodo-downloaded data folder
data_folder <- "~/Downloads/sushie_real_data_results/all_results/"


rnaseq_cov <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE")

rnaseq_indep <- read_tsv(glue("{data_folder}/rnaseq_indep.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE-Indep")

rnaseq_meta <- read_tsv(glue("{data_folder}/rnaseq_normal.meta_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "Meta-SuSiE")

rnaseq_mega <- read_tsv(glue("{data_folder}/rnaseq_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuSiE")

rnaseq_mesusie <- read_tsv(glue("{data_folder}/rnaseq_mesusie_cs.tsv.gz")) %>%
  distinct(trait) %>%
  mutate(name = "MESuSiE")

proteins_cov <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE")

proteins_indep <- read_tsv(glue("{data_folder}/proteins_indep.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE-Indep")

proteins_meta <- read_tsv(glue("{data_folder}/proteins_normal.meta_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "Meta-SuSiE")

proteins_mega <- read_tsv(glue("{data_folder}/proteins_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuSiE")

proteins_mesusie <- read_tsv(glue("{data_folder}/proteins_mesusie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "MESuSiE")

genoa_cov <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE")

genoa_indep <- read_tsv(glue("{data_folder}/genoa_indep.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuShiE-Indep")

genoa_meta <- read_tsv(glue("{data_folder}/genoa_normal.meta_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "Meta-SuSiE")

genoa_mega <- read_tsv(glue("{data_folder}/genoa_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "SuSiE")

genoa_mesusie <- read_tsv(glue("{data_folder}/genoa_mesusie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(name = "MESuSiE")

rnaseq_enet <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "Elastic Net") %>%
  distinct(trait, name) %>%
  filter(trait %in% rnaseq_her$trait)

rnaseq_lasso <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "LASSO") %>%
  distinct(trait, name) %>%
  filter(trait %in% rnaseq_her$trait)

rnaseq_gblup <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "gBLUP") %>%
  distinct(trait, name) %>%
  filter(trait %in% rnaseq_her$trait)

proteins_enet <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "Elastic Net") %>%
  distinct(trait, name) %>%
  filter(trait %in% proteins_her$trait)

proteins_lasso <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "LASSO") %>%
  distinct(trait, name) %>%
  filter(trait %in% proteins_her$trait)

proteins_gblup <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "gBLUP") %>%
  distinct(trait, name) %>%
  filter(trait %in% proteins_her$trait)

genoa_enet <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "Elastic Net") %>%
  distinct(trait, name) %>%
  filter(trait %in% genoa_her$trait)

genoa_lasso <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "LASSO") %>%
  distinct(trait, name) %>%
  filter(trait %in% genoa_her$trait)

genoa_gblup <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz")) %>%
  mutate(name = "gBLUP") %>%
  distinct(trait, name) %>%
  filter(trait %in% genoa_her$trait)

df_all <- bind_rows(
  rnaseq_cov, rnaseq_indep, rnaseq_meta, rnaseq_mega, rnaseq_mesusie,
  proteins_cov, proteins_indep, proteins_meta, proteins_mega, proteins_mesusie,
  genoa_cov, genoa_indep, genoa_meta, genoa_mega, genoa_mesusie,
  rnaseq_enet, rnaseq_lasso, rnaseq_gblup,
  proteins_enet, proteins_lasso, proteins_gblup,
  genoa_enet, genoa_lasso, genoa_gblup)

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

twas_colors <- c(method_colors,
  "LASSO" = "#fb8072", "Elastic Net" = "#b3de69", "gBLUP" = "#fccde5")

rnaseq_r2 <- read_tsv(glue("{data_folder}/rnaseq_sushie.r2.tsv.gz")) %>%
  mutate(study = "mesa.mrna") %>%
  pivot_longer(sushie:cross) %>%
  bind_rows(
    read_tsv(glue("{data_folder}/rnaseq_mesusie.r2.tsv.gz")) %>%
      mutate(study = "mesa.mrna",
        name = "mesusie",
        type = "r2") %>%
      select(type, trait, study, name, value = mesusie)
  )

proteins_r2 <- read_tsv(glue("{data_folder}/proteins_sushie.r2.tsv.gz")) %>%
  mutate(study = "mesa.proteins") %>%
  pivot_longer(sushie:cross) %>%
  bind_rows(
    read_tsv(glue("{data_folder}/proteins_mesusie.r2.tsv.gz")) %>%
      mutate(study = "mesa.proteins",
        name = "mesusie",
        type = "r2") %>%
      select(type, trait, study, name, value = mesusie)
  )

genoa_r2 <- read_tsv(glue("{data_folder}/genoa_sushie.r2.tsv.gz")) %>%
  mutate(study = "genoa.mrna") %>%
  pivot_longer(sushie:cross) %>%
  bind_rows(
    read_tsv(glue("{data_folder}/genoa_mesusie.r2.tsv.gz")) %>%
      mutate(study = "genoa.mrna",
        name = "mesusie",
        type = "r2") %>%
      select(type, trait, study, name, value = mesusie)
  )

rnaseq_corr <-
  read_tsv(glue("{data_folder}/rnaseq_corr.tsv.gz"))

proteins_corr <-
  read_tsv(glue("{data_folder}/proteins_corr.tsv.gz"))

genoa_corr <- read_tsv(glue("{data_folder}/genoa_corr.tsv.gz"))

df_r2 <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(name != "cross") %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "mesusie", "enet",
      "lasso",  "gblup"),
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE",
      "MESuSiE", "Elastic Net", "LASSO", "gBLUP")),
    study = factor(study,
      levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Proteins",
        "GENOA mRNA"))) %>%
  inner_join(df_all)

df_r2 %>%
  group_by(name) %>%
  summarize(mval = mean(value),
    se = sd(value)/sqrt(n())) %>%
  arrange(desc(mval))

r2_res <- tibble()
for (met in c("SuShiE-Indep",  "Meta-SuSiE", "SuSiE",  "MESuSiE",
  "Elastic Net", "LASSO", "gBLUP")) {
  df_tmp <- df_r2 %>%
    filter(name %in% c("SuShiE", met)) %>%
    group_by(trait, study) %>%
    filter(n() == 2) %>%
    mutate(name = factor(name, levels = c("SuShiE", met)))
  
  r2_res <- r2_res %>%
    bind_rows(tidy(lm(value ~ name + study, df_tmp))[2,])
}

r2_res
r2_res %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

r2_res_plot <- tibble()
for (met in c("SuShiE-Indep",  "Meta-SuSiE", "SuSiE",  "MESuSiE",
  "Elastic Net", "LASSO", "gBLUP")) {
  r2_res_plot <- r2_res_plot %>%
    bind_rows(
      df_r2 %>%
    filter(name %in% c("SuShiE", met)) %>%
    group_by(trait, study) %>%
    filter(n() == 2) %>%
    mutate(name = factor(name, levels = c("SuShiE", met),
      labels = c("SuShiE", "Other Method")),
      type = met)
    )
}

rnaseq_cov_original <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz"))

proteins_cov_original <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz"))

genoa_cov_original <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz"))

heter_genes <- bind_rows(
  rnaseq_corr %>%
    inner_join(rnaseq_cov_original %>%
        filter(!is.na(snp)) %>%
        distinct(trait, CSIndex)) %>%
    select(trait, CSIndex, contains("corr")) %>%
    pivot_longer(cols = 3:5) %>%
    rename(type = name, corr = value) %>%
    mutate(study = "mesa.mrna") %>%
    select(study, trait, type, CSIndex, corr) %>%
    filter(corr < 0.9) %>%
    distinct(trait),
  proteins_corr %>%
    inner_join(proteins_cov_original %>%
        filter(!is.na(snp)) %>%
        distinct(trait, CSIndex)) %>%
    select(trait, CSIndex, contains("corr")) %>%
    pivot_longer(cols = 3:5) %>%
    rename(type = name, corr = value) %>%
    mutate(study = "mesa.mrna") %>%
    select(study, trait, type, CSIndex, corr) %>%
    filter(corr < 0.9) %>%
    distinct(trait),
  genoa_corr %>%
    inner_join(genoa_cov_original %>%
        filter(!is.na(snp)) %>%
        distinct(trait, CSIndex)) %>%
    select(trait, CSIndex, contains("corr")) %>%
    pivot_longer(cols = 3) %>%
    rename(type = name, corr = value) %>%
    mutate(study = "mesa.mrna") %>%
    select(study, trait, type, CSIndex, corr) %>%
    filter(corr < 0.9) %>%
    distinct(trait)
)

df_r2_heter <- df_r2 %>%
  filter(trait %in% heter_genes$trait)

r2_res_heter <- tibble()
for (met in c("SuShiE-Indep",  "Meta-SuSiE", "SuSiE",  "MESuSiE",
  "Elastic Net", "LASSO", "gBLUP")) {
  df_tmp <- df_r2_heter %>%
    filter(name %in% c("SuShiE", met)) %>%
    group_by(trait, study) %>%
    filter(n() == 2) %>%
    mutate(name = factor(name, levels = c("SuShiE", met)))
  
  r2_res_heter <- r2_res_heter %>%
    bind_rows(tidy(lm(value ~ name + study, df_tmp))[2,])
}

r2_res_heter %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

r2_res_heter_plot <- r2_res_plot %>%
  filter(trait %in% heter_genes$trait)

df_r2_comp <- r2_res_plot %>%
  group_by(name, type) %>%
  summarize(mval = mean(value),
    se_upp = mval + 1.96 * sd(value)/sqrt(n()),
    se_low = mval - 1.96 * sd(value)/sqrt(n())) %>%
  mutate(class = "all") %>%
  bind_rows(r2_res_heter_plot %>%
      group_by(name, type) %>%
      summarize(mval = mean(value),
        se_upp = mval + 1.96 * sd(value)/sqrt(n()),
        se_low = mval - 1.96 * sd(value)/sqrt(n())) %>%
      mutate(class = "heter")) %>%
  mutate(class = factor(class, levels = c("all", "heter"),
    labels = c("A: All e/pGenes", "B:e/pGenes Exhibited Heterogeneity")),
    type = factor(type , levels = c("SuShiE-Indep", "Meta-SuSiE",
      "SuSiE", "MESuSiE", "Elastic Net", "LASSO", "gBLUP")))

r2_colors <-c("SuShiE" = "#1b9e77", "Other Method" = "#8856a7")
ggplot(df_r2_comp,
  aes(x = type, y = mval, color = name)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = se_low, ymax = se_upp),
    position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(cols = vars(class), scales="free_y") +
  scale_color_manual(values = r2_colors) +
  ylab("Average CV r-squared") +
  xlab("Prediction Method") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 8, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(size = 8, face="bold", angle=45, vjust = 0.5),
    title = element_text(size = 10, face="bold"),
    axis.text.y =element_text(size = 8, face="bold"))

# ggsave("./plots/s25.png", width = p_width, height = p_height+2)

df_r2_cross <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(name %in% c("sushie", "cross")) %>%
  group_by(name) %>%
  mutate(name = factor(name,
    levels = c("sushie", "cross"),
    labels = c("Ancestry-matched weights", "Cross-ancestry weights"))) %>%
  group_by(name) %>%
  summarize(mval = mean(value),
    se_upp = mval + 1.96 * sd(value)/sqrt(n()),
    se_low = mval - 1.96 * sd(value)/sqrt(n())) 

cross_colors <- c("#1b9e77", "#ff7f00")

ggplot(df_r2_cross, aes(x=name, y = mval, color=name)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = se_low, ymax = se_upp),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_color_manual(values = cross_colors) +
  ylab("Average CV r-squared") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 8, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    axis.title.x=element_blank(),
    title = element_text(size = 10, face="bold"),
    axis.text=element_text(size = 8, face="bold"),
    text=element_text(size = 8))

# ggsave("./plots/s26.png", width = p_width/2, height = p_height)

