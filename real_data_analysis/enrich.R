# this script is to create main figure 3
library(tidyverse)
library(broom)
library(RColorBrewer)
library(glue)
library(ggpubr)
source("./utils.R")

text_size <- 4
title_size <- 9
legend_size <- 8
bar_width <- 0.5

# rnaseq
rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.sushie_cs.tsv.gz")

proteins_cov<- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.sushie_cs.tsv.gz")

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.sushie_cs.tsv.gz")


# qtl
rnaseq_qtl <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(new_n = ifelse(n < 6, n, "6 or more")) %>%
  group_by(new_n) %>%
  mutate(ct = n()) %>%
  ungroup() %>%
  mutate(perc_val = round((ct/n() * 100), 2),
    perc_str = paste0(perc_val,"%")) %>%
  distinct(n,new_n, ct, perc_str,perc_val) %>%
  mutate(study = "TOPMed-MESA mRNA")  %>%
  filter(n <= 6)

proteins_qtl <- proteins_cov%>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(new_n = ifelse(n < 6, n, "6 or more")) %>%
  group_by(new_n) %>%
  mutate(ct = n()) %>%
  ungroup() %>%
  mutate(perc_val = round((ct/n() * 100), 2),
    perc_str = paste0(perc_val,"%")) %>%
  distinct(n,new_n, ct, perc_str,perc_val) %>%
  mutate(study = "TOPMed-MESA Protein")  %>%
  filter(n <= 6)

genoa_qtl <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  summarize(n = length(unique(CSIndex))) %>%
  mutate(new_n = ifelse(n < 6, n, "6 or more")) %>%
  group_by(new_n) %>%
  mutate(ct = n()) %>%
  ungroup() %>%
  mutate(perc_val = round((ct/n() * 100), 2),
    perc_str = paste0(perc_val,"%")) %>%
  distinct(n,new_n, ct, perc_str,perc_val) %>%
  mutate(study = "GENOA mRNA")  %>%
  filter(n <= 6)

# df_qtl <- bind_rows(rnaseq_qtl, genoa_qtl, proteins_qtl) %>%
df_qtl <- rnaseq_qtl %>%
  bind_rows(proteins_qtl, genoa_qtl) %>%
  mutate(new_n = factor(new_n, levels = c(1:5, "6 or more"),
    labels = paste0(c(1:5, "6 or more"), "\n")),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

main_theme_p1 <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    # legend.position = "none",
    strip.text = element_blank(), # Removes facet labels
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 8),
    axis.text = element_text(face = "bold", size = 9))
}

mp1 <- ggplot(df_qtl, aes(x=new_n, y=perc_val, fill = study)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = main_study_color) +
  # geom_text(aes(label = perc_str), vjust = -0.3, size = 1.5) +
  xlab("Number of molQTLs") +
  ylab("Gene Percentage") +
  scale_y_continuous(breaks=seq(0, 60, 10), limits = c(0, 65)) +
  main_theme_p1()

leg <- get_legend(mp1, position = "bottom")

rnaseq_tss <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_tss.tsv.gz") %>%
  group_by(Type) %>%
  mutate(n_gene = length(unique(gene))) %>%
  group_by(bins, Type) %>%
  summarize(avg_signal = mean(mean),
    n = n(),
    n_gene = mean(n_gene)) %>%
  mutate(study = "TOPMed-MESA mRNA")

proteins_tss <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_tss.tsv.gz") %>%
  group_by(Type) %>%
  mutate(n_gene = length(unique(gene))) %>%
  group_by(bins, Type) %>%
  summarize(avg_signal = mean(mean),
    n = n(),
    n_gene = mean(n_gene)) %>%
  mutate(study = "TOPMed-MESA Protein")

genoa_tss <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_tss.tsv.gz") %>%
  group_by(Type) %>%
  mutate(n_gene = length(unique(gene))) %>%
  group_by(bins, Type) %>%
  summarize(avg_signal = mean(mean),
    n = n(),
    n_gene = mean(n_gene)) %>%
  mutate(study = "GENOA mRNA")

df_pip <- rnaseq_tss %>%
  filter(Type %in% "pip_all") %>%
  bind_rows(proteins_tss %>%
      filter(Type %in% "pip_all")) %>%
  bind_rows(genoa_tss %>%
      filter(Type %in% "pip_all")) %>%
  mutate(bins = ((bins - 1000)*500)/1000,
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

mp2 <- ggplot(df_pip, aes(x = bins, y = avg_signal, color=study)) +
  geom_line() +
  scale_color_manual(values = main_study_color) +
  facet_grid(rows = vars(study)) +
  xlab("Distance to TSS (kb)") +
  ylab("Average PIP") +
  scale_y_continuous(breaks=c(0, 0.03, 0.06), limits = c(0,0.06)) +
  scale_x_continuous(breaks = c(-500, -250, 0, 250, 500),
    labels = c("-500\n", "-250\n", "0\n", "250\n", "500\n")) + 
  main_theme_p1()

# enrich
rnaseq_enrich <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_enrich_all.tsv.gz") %>%
  filter(!anno %in% c("tss_all", "tss_protein"))

proteins_enrich <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_enrich_all.tsv.gz") %>%
  filter(!anno %in% c("tss_all", "tss_protein"))

genoa_enrich <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_enrich_all.tsv.gz") %>%
  filter(!anno %in% c("tss_all", "tss_protein"))

enrich_list <- unique(rnaseq_enrich$anno)
ldsc_list <- enrich_list[grepl("LDSC", enrich_list)]
ldsc_list <- ldsc_list[order(ldsc_list)]
encode_list <- c("PLS", "pELS", "dELS", "CTCF-bound", "DNase-H3K4me3")
atac_list <- enrich_list[!enrich_list %in% c(ldsc_list,
  encode_list, "snATAC-seq-frozen-peaks", "megakaryocyte")]

rnaseq_anno <- rnaseq_enrich %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% "pip_cov") %>%
  filter(anno %in% encode_list) %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  mutate(study = "TOPMed-MESA mRNA")

proteins_anno <- proteins_enrich %>%
  filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% "pip_cov") %>%
  filter(anno %in% encode_list) %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  mutate(study = "TOPMed-MESA Protein")

genoa_anno <- genoa_enrich %>%
  filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% "pip_cov") %>%
  filter(anno %in% encode_list) %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  mutate(study = "GENOA mRNA")

df_anno <- rnaseq_anno %>%
  bind_rows(proteins_anno, genoa_anno) %>%
  select(anno, log_est, low_int, upp_int, study) %>%
  mutate(anno = factor(anno,
    levels = c("PLS", "pELS", "dELS",  "CTCF-bound", "DNase-H3K4me3"),
    # labels = c("Promoter", "Proximal\nEnhancer",
    #   "Distal\nEnhancer", "CTCF", "DNase\nH3K4me3")),
    labels = c("PLS", "pELS",
      "dELS", "CTCF", "DNase\nH3K4me3")),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

mp3 <- ggplot(df_anno,
  aes(x = reorder(anno, log_est, decreasing = TRUE),
    y = log_est, color=study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_int, ymax = upp_int),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab(bquote(bold(Log[2](Enrichment)))) +
  xlab("Functional Annotations") +
  scale_y_continuous(limits = c(-0.6,1.05), breaks = c(-0.5, 0, 0.5, 1)) +
  scale_color_manual(values = main_study_color) +
  geom_hline(yintercept=0, linetype="dashed") +
  main_theme_p1()


enrich_res <- rnaseq_enrich %>%
  filter(anno %in% c(atac_list, "snATAC-seq-frozen-peaks")) %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  mutate(study = "TOPMed-MESA mRNA") %>%
  bind_rows(proteins_enrich %>%
      filter(anno %in% c(atac_list, "snATAC-seq-frozen-peaks")) %>%
      filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov") %>%
      mutate(study = "TOPMed-MESA Protein"),
    genoa_enrich %>%
      filter(anno %in% c(atac_list, "snATAC-seq-frozen-peaks")) %>%
      filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov") %>%
      mutate(study = "GENOA mRNA")) %>%
  group_by(anno, study) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight) / sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96 * meta_se) / log(2),
    upp_int = (meta_est + 1.96 * meta_se) / log(2)) %>%
  mutate(anno = factor(anno,
    levels = c("regulatory-T",
      "memory-CD-T", "adaptive-NK", "activated-CD-T",
      "naive-T", "naive-B", "memory-B",
      "non-classical-monocyte",
      "classical-monocyte", "snATAC-seq-frozen-peaks"),
    labels = c("Regulatory\nT", "Memeory\nCD-T",
      "Adaptive\nNK",  "Activated\nCD-T", "Naive\nT", "Naive\nB", "Memory\nB",
      "Non-classical\nMonocyte",  "Classical\nMonocyte", "PBMC")),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

main1 <- ggarrange(mp1, mp2, mp3, align = "h", nrow=1,
  labels = c("A", "B", "C"), font.label = list(size = 10),
  legend.grob =leg, legend = "bottom")

# ggsave("./plots/p3.png", width = p_width, height = p_height+1)

mp4 <- ggplot(enrich_res,
  aes(x = reorder(anno, log_est, decreasing = TRUE),
    y = log_est, color=study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_int, ymax = upp_int),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab(bquote(bold(Log[2](Enrichment)))) +
  xlab("Cell-type/tissue-specific cCREs by sn/scATAC-seq") +
  scale_y_continuous(limits = c(-0.2, 0.8), breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8)) +
  scale_color_manual(values = main_study_color) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_blank(), # Removes facet labels
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title.y = element_text(face="bold", size = 8),
    axis.title.x = element_text(face="bold", size = 8,
      margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text = element_text(face = "bold", size = 8))

# ggsave("./plots/s12.png", , width = p_width-1, height = p_height)

# alpha tss
df_pip <- rnaseq_tss %>%
  filter(Type %in% paste0("alpha_l", 1:6)) %>%
  bind_rows(proteins_tss %>%
      filter(Type %in% paste0("alpha_l", 1:6))) %>%
  bind_rows(genoa_tss %>%
      filter(Type %in% paste0("alpha_l", 1:6))) %>%
  mutate(bins = ((bins - 1000)*500)/1000,
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")),
    Type = factor(Type, levels = paste0("alpha_l", 1:6),
      labels = paste0("L", 1:6)))

ggplot(df_pip, aes(x = bins, y = avg_signal, color=study)) +
  geom_line() +
  scale_color_manual(values = main_study_color) +
  facet_grid(cols = vars(study), rows = vars(Type), scales = "free",) +
  # facet_wrap(study ~ Type, scales = "free", ncol=3, dir = "v") +
  xlab("Distance to TSS (kb)") +
  ylab("Average Posterior Probability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text.x = element_blank(), # Removes facet labels
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 8),
    axis.text = element_text(face = "bold", size = 8))

# ggsave("./plots/s13.png", , width = p_width-1.5, height = p_height+2)

# single effect enrichment
rnaseq_alpha <- rnaseq_enrich %>%
  filter(converged == 1 & reg == 0) %>%
  filter(!grepl("pip", method) & !grepl("cs", method)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  filter(anno %in% c(encode_list, atac_list, "snATAC-seq-frozen-peaks")) %>%
  mutate(study = "TOPMed-MESA mRNA")

proteins_alpha <- proteins_enrich %>%
  filter(converged == 1 & reg == 0) %>%
  filter(!grepl("pip", method) & !grepl("cs", method)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  filter(anno %in% c(encode_list, atac_list, "snATAC-seq-frozen-peaks")) %>%
  mutate(study = "TOPMed-MESA Protein")

genoa_alpha <- rnaseq_enrich %>%
  filter(converged == 1 & reg == 0) %>%
  filter(!grepl("pip", method) & !grepl("cs", method)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_est = sum(est*weight)/sum(weight),
    meta_se = 1 / sqrt(sum(weight)),
    log_est = meta_est / log(2),
    low_int = (meta_est - 1.96*meta_se)/log(2),
    upp_int = (meta_est + 1.96*meta_se)/log(2),
    n = n()) %>%
  filter(anno %in% c(encode_list, atac_list, "snATAC-seq-frozen-peaks")) %>%
  mutate(study = "GENOA mRNA")

anno_alpha_encode <- rnaseq_alpha %>%
  bind_rows(proteins_alpha, genoa_alpha) %>%
  filter(method %in% paste0("cov_l", 1:6)) %>%
  mutate(log_est = meta_est / log(2)) %>%
  filter(anno %in% c("PLS", "pELS", "dELS",  "CTCF-bound", "DNase-H3K4me3")) %>%
  mutate(anno = factor(anno,
    levels = c("PLS", "pELS", "dELS",  "CTCF-bound", "DNase-H3K4me3"),
    labels = c("Promoter", "Proximal\nEnhancer",
      "Distal\nEnhancer", "CTCF", "DNase\nH3K4me3")),
    method = factor(method, levels = paste0("cov_l", 1:6),
      labels = paste0("L", 1:6)),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

anno_alpha_chiou <- rnaseq_alpha %>%
  bind_rows(proteins_alpha, genoa_alpha) %>%
  filter(method %in% paste0("cov_l", 1:6)) %>%
  mutate(log_est = meta_est / log(2)) %>%
  filter(anno %in% c("regulatory-T",
    "memory-CD-T", "adaptive-NK", "activated-CD-T",
    "naive-T", "naive-B", "memory-B",
    "non-classical-monocyte",
    "classical-monocyte", "snATAC-seq-frozen-peaks")) %>%
  mutate(anno = factor(anno,
    levels = c("regulatory-T",
      "memory-CD-T", "adaptive-NK", "activated-CD-T",
      "naive-T", "naive-B", "memory-B",
      "non-classical-monocyte",
      "classical-monocyte", "snATAC-seq-frozen-peaks"),
    labels = c("Regulatory\nT", "Memeory\nCD-T",
      "Adaptive\nNK",  "Activated\nCD-T", "Naive\nT", "Naive\nB", "Memory\nB",
      "Non-classical\nMonocyte",  "Classical\nMonocyte", "PBMC")),
    method = factor(method, levels = paste0("cov_l", 1:6),
      labels = paste0("L", 1:6)),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")))

m1 <- ggplot(anno_alpha_encode,
  aes(x = reorder(anno, log_est, decreasing = TRUE),
    y = log_est, color=study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_int, ymax = upp_int),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab(bquote(bold(Log[2](Enrichment)))) +
  xlab("cCREs by ATAC-seq in ENCODE") +
  facet_grid(rows = vars(method), scales = "free", switch = "both") +
  # scale_y_continuous(limits = c(-0.5,1.05), breaks = c(-0.5, 0, 0.5, 1)) +
  scale_color_manual(values = main_study_color) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text = element_text(face="bold", size = 8), 
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 8),
    axis.text = element_text(face = "bold", size = 8))

m2 <- ggplot(anno_alpha_chiou,
  aes(x = reorder(anno, log_est, decreasing = TRUE),
    y = log_est, color=study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_int, ymax = upp_int),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab(bquote(bold(Log[2](Enrichment)))) +
  xlab("cCREs by sn/sc ATAC-seq in Chiou et al. and Satpathy et al.") +
  facet_grid(rows = vars(method), scales = "free", switch = "both") +
  # scale_y_continuous(limits = c(-0.5,1.05), breaks = c(-0.5, 0, 0.5, 1)) +
  scale_color_manual(values = main_study_color) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title.x=element_text(face="bold", size = 8),
    axis.title.y=element_blank(),
    axis.text = element_text(face = "bold", size = 8))

ggarrange(m1, m2, nrow=1, widths=c(1,2), align = "h", common.legend = TRUE,
  legend = "bottom")

# ggsave("./plots/s14.png", , width = p_width+3, height = p_height+5)


# creating tables

rnaseq_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_indep.sushie_cs.tsv.gz")

rnaseq_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.meta_cs.tsv.gz")

rnaseq_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.mega_cs.tsv.gz")

rnaseq_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_susiex_cs.tsv.gz")

rnaseq_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_mesusie_cs.tsv.gz")

proteins_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_indep.sushie_cs.tsv.gz")

proteins_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.meta_cs.tsv.gz")

proteins_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_normal.mega_cs.tsv.gz")

proteins_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_susiex_cs.tsv.gz")

proteins_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/proteins_mesusie_cs.tsv.gz")

genoa_indep <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_indep.sushie_cs.tsv.gz")

genoa_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.meta_cs.tsv.gz")

genoa_mega <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_normal.mega_cs.tsv.gz")

genoa_susiex <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_susiex_cs.tsv.gz")

genoa_mesusie <- read_tsv("~/Documents/github/data/sushie_results/real2/genoa_mesusie_cs.tsv.gz")


meta_anno <- tibble(Study = "LDSC", Annotation = ldsc_list) %>%
  bind_rows(tibble(Study = "ENCODE", Annotation = encode_list),
    tibble(Study = "Chiou et al. (snATAC-seq)", Annotation = atac_list),
    tibble(Study = "Satpathy et al. (scATAC-seq)",
      Annotation = "snATAC-seq-frozen-peaks"))

# aggregate for manuscript across studies
enrich_res <- rnaseq_enrich %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  bind_rows(proteins_enrich %>%
      filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov"),
    genoa_enrich %>%
      filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov")) %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_all <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

# output rnaseq
enrich_res <- rnaseq_enrich %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_rnaseq <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

# output protiens
enrich_res <- proteins_enrich %>%
  filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_proteins <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

# output genoa
enrich_res <- genoa_enrich %>%
  filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_genoa <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

table_enrich <- table_enrich_all %>%
  left_join(table_enrich_rnaseq, by = c("Study", "Annotation")) %>%
  left_join(table_enrich_proteins, by = c("Study", "Annotation")) %>%
  left_join(table_enrich_genoa, by = c("Study", "Annotation"))

# write_tsv(table_enrich %>% select(-Study, - Annotation), "./tables/s5.tsv", col_names=FALSE)

# in-text number 
enrich_res <- rnaseq_enrich %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method == "pip_cov") %>%
  bind_rows(proteins_enrich %>%
      filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov"),
    genoa_enrich %>%
      filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
      filter(converged == 1 & reg == 0) %>%
      filter(method == "pip_cov")) %>%
  filter(anno %in% c(ldsc_list, atac_list, encode_list,
    "snATAC-seq-frozen-peaks")) %>%
  group_by(anno) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

enrich_res %>%
  filter(meta_z > qnorm(0.05/89, lower.tail = FALSE))

enrich_res %>%
  filter(anno %in% encode_list)

enrich_res %>%
  filter(anno %in% c(atac_list, "snATAC-seq-frozen-peaks"))

# method comparison
comp_enrich <- function(df, method_name) {
  tmp_df <- df %>%
    filter(converged == 1 & reg == 0) %>%
    filter(anno %in% c(ldsc_list, atac_list, encode_list,
      "snATAC-seq-frozen-peaks")) %>%
    filter(method %in% c("pip_cov", method_name)) %>%
    group_by(trait, anno) %>%
    filter(n() == 2) %>%
    group_by(method) %>%
    mutate(weight = 1 / (se*2)) %>%
    summarize(est = sum(est*weight)/sum(weight),
      se = 1 / sqrt(sum(weight))) %>%
    mutate(method = ifelse(method == "pip_cov", "SuShiE", "Comp"),
      type = method_name) %>%
    pivot_wider(names_from = method, values_from = c(est, se)) %>%
    mutate(compz = (est_SuShiE - est_Comp) /
        sqrt(se_SuShiE^2 + se_Comp^2)) %>%
    select(type, est_SuShiE, se_SuShiE, est_Comp, se_Comp, compz)
  
  return(tmp_df)
}


enrich_all <- tibble()
for (method in c("pip_indep", "pip_meta", "susiex_cs", "mesusie_cs",
  "pip_mega")) {
  enrich_all <- enrich_all %>%
    bind_rows(
      comp_enrich(bind_rows(rnaseq_enrich,
        proteins_enrich, genoa_enrich), method) %>%
        mutate(dataset = "across"),
      comp_enrich(rnaseq_enrich, method) %>%
        mutate(dataset = "TOPMed-MESA mRNA"),
      comp_enrich(proteins_enrich, method) %>%
        mutate(dataset = "TOPMed-MESA Protein"),
      comp_enrich(genoa_enrich, method) %>%
        mutate(dataset = "GENOA mRNA")
    )
}


df_enrich_comp <- enrich_all %>%
  mutate(dataset = factor(dataset,levels = c("across",
    "TOPMed-MESA mRNA", "TOPMed-MESA Protein", "GENOA mRNA")),
    type = factor(type, levels = c("pip_indep",
      "pip_meta", "pip_mega", "susiex_cs", "mesusie_cs"),
      labels = c("SuShiE vs SuShiE-Indep",
        "SuShiE vs Meta-SuSiE",
        "SuShiE vs SuSiE",
        "SuShiE vs SuSiEx",
        "SuShiE vs MESuSiE"))) %>%
  select(dataset, type, est_SuShiE, se_SuShiE, est_Comp, se_Comp, compz) %>%
  arrange(dataset, type)

# write_tsv(df_enrich_comp %>% select(-dataset, -type), "./tables/s6.tsv", col_names=FALSE)


# single effect alpha
# aggregate for manuscript across studies
enrich_res <- rnaseq_enrich %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% paste0("cov_l", 1:10)) %>%
  bind_rows(proteins_enrich %>%
      filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
      filter(method %in% paste0("cov_l", 1:10)) %>%
      filter(method == "pip_cov"),
    genoa_enrich %>%
      filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
      filter(method %in% paste0("cov_l", 1:10)) %>%
      filter(method == "pip_cov")) %>%
  filter(anno %in% c(ldsc_list, atac_list, encode_list,
    "snATAC-seq-frozen-peaks")) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_all <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_z))

# output rnaseq
enrich_res <- rnaseq_enrich %>%
  filter(anno %in% c(ldsc_list, atac_list, encode_list,
    "snATAC-seq-frozen-peaks")) %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% paste0("cov_l", 1:10)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_rnaseq <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

# output protiens
enrich_res <- proteins_enrich %>%
  filter(anno %in% c(ldsc_list, atac_list, encode_list,
    "snATAC-seq-frozen-peaks")) %>%
  filter(trait %in% filter(proteins_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% paste0("cov_l", 1:10)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_proteins <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

# output genoa
enrich_res <- genoa_enrich %>%
  filter(anno %in% c(ldsc_list, atac_list, encode_list,
    "snATAC-seq-frozen-peaks")) %>%
  filter(trait %in% filter(genoa_cov, !is.na(snp))$trait) %>%
  filter(converged == 1 & reg == 0) %>%
  filter(method %in% paste0("cov_l", 1:10)) %>%
  group_by(anno, method) %>%
  mutate(weight = 1 / (se*2)) %>%
  summarize(meta_odds = exp(sum(est*weight)/sum(weight)),
    meta_z = sum(est*weight) / sqrt(sum(weight)),
    n = n())

table_enrich_genoa <- meta_anno %>%
  left_join(enrich_res %>%
      rename(Annotation = anno),
    by = "Annotation") %>%
  mutate(Study = factor(Study, levels = c("ENCODE", "Chiou et al. (snATAC-seq)",
    "Satpathy et al. (scATAC-seq)", "LDSC")),
    Annotation = ifelse(Annotation == "snATAC-seq-frozen-peaks", "PBMC",
      ifelse(Study == "LDSC", gsub("LDSC_", "", Annotation), Annotation))) %>%
  arrange(Study, desc(meta_odds))

table_enrich <- table_enrich_all %>%
  left_join(table_enrich_rnaseq, by = c("Study", "Annotation", "method")) %>%
  left_join(table_enrich_proteins, by = c("Study", "Annotation", "method")) %>%
  left_join(table_enrich_genoa, by = c("Study", "Annotation", "method")) %>%
  mutate(method = gsub("cov_l", "L", method))

# write_tsv(table_enrich, "./tables/s7.tsv", col_names=FALSE)

rnaseq_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_rnaseq_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X12, TSS = X6)

proteins_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_proteins_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X15, TSS = X6)

genoa_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/genoa_sushie_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X2, TSS = X6)

rnaseq_dist <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, CSIndex, pip_all, trait) %>%
  left_join(rnaseq_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait, CSIndex) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  group_by(CSIndex) %>%
  summarize(mval = mean(dist),
    se = sd(dist)/sqrt(n()),
    se_upp = mval+se*1.96,
    se_low = mval-se*1.96) %>%
  mutate(study = "mesa.mrna")

proteins_dist <- proteins_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, CSIndex, pip_all, trait) %>%
  left_join(proteins_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait, CSIndex) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  group_by(CSIndex) %>%
  summarize(mval = mean(dist),
    se = sd(dist)/sqrt(n()),
    se_upp = mval+se*1.96,
    se_low = mval-se*1.96) %>%
  mutate(study = "mesa.proteins")

genoa_dist <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, CSIndex, pip_all, trait) %>%
  left_join(genoa_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait, CSIndex) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  group_by(CSIndex) %>%
  summarize(mval = mean(dist),
    se = sd(dist)/sqrt(n()),
    se_upp = mval+se*1.96,
    se_low = mval-se*1.96) %>%
  mutate(study = "genoa.mrna")

total_dist <- bind_rows(rnaseq_dist,
  proteins_dist, genoa_dist) %>%
  mutate(study = factor(study,
    levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
    labels = c("TOPMed-MESA mRNA", "TOPMed-MESA Proteins",
      "GENOA mRNA")),
    CSIndex = factor(CSIndex, levels = 1:10))


ggplot(total_dist, aes(x = CSIndex, y = mval, color = study)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = se_low, ymax = se_upp),
    position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(rows=vars(study)) +
  scale_color_manual(values = main_study_color) +
  xlab("cis-molQTL Index") +
  ylab("Expected Distance to TSS") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text = element_text(face="bold", size = 8), 
    panel.spacing = unit(0.5, "lines"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 8),
    axis.text = element_text(face = "bold", size = 8))

# ggsave("./plots/s15.png", width = p_width-3, height = p_height+3)


two_dist <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, CSIndex, pip_all, trait) %>%
  left_join(rnaseq_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait, CSIndex) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  bind_rows(
    proteins_dist <- proteins_cov %>%
      filter(!is.na(snp)) %>%
      select(snp, pos, CSIndex, pip_all, trait) %>%
      left_join(proteins_ref, by = "trait") %>%
      mutate(dist = abs(pos - TSS) + 1) %>%
      group_by(trait, CSIndex) %>%
      summarize(dist = sum(pip_all*dist) / sum(pip_all)),
    genoa_dist <- genoa_cov %>%
      filter(!is.na(snp)) %>%
      select(snp, pos, CSIndex, pip_all, trait) %>%
      left_join(genoa_ref, by = "trait") %>%
      mutate(dist = abs(pos - TSS) + 1) %>%
      group_by(trait, CSIndex) %>%
      summarize(dist = sum(pip_all*dist) / sum(pip_all))
  ) %>%
  mutate(type = ifelse(CSIndex %in% 1:3, "first", "last")) %>%
  group_by(type) %>%
  summarize(mval = mean(dist),
    se = sd(dist)/sqrt(n()),
    se_upp = mval+se*1.96,
    se_low = mval-se*1.96)

(two_dist$mval[1] - two_dist$mval[2])/sqrt(two_dist$se[1]^2 + two_dist$se[2]^2)

pnorm(abs((two_dist$mval[1] - two_dist$mval[2])
  /sqrt(two_dist$se[1]^2 + two_dist$se[2]^2)), lower.tail = FALSE)
