library(tidyverse)
library(broom)
library(ggpubr)

source("./utils.R")

rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz") 

proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz") 

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz")

rnaseq_her <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_her.tsv.gz")

proteins_her <-
  read_tsv("~/Documents/github/data/sushie_results/real/proteins_her.tsv.gz")

genoa_her <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_her.tsv.gz")

# get legend
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
    axis.text = element_text(face = "bold", size = 7))
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


rnaseq_genes <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait) 

proteins_genes <- proteins_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait) 

genoa_genes <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait) 

rnaseq_sig_trait <- rnaseq_her %>%
  filter(trait %in% rnaseq_genes$trait) %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==3) %>%
  distinct(trait)

proteins_sig_trait <- proteins_her %>%
  filter(trait %in% proteins_genes$trait) %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==3) %>%
  distinct(trait)

genoa_sig_trait <- genoa_her %>%
  filter(trait %in% genoa_genes$trait) %>%
  filter(p_value < 0.05) %>%
  group_by(trait) %>%
  filter(n()==2) %>%
  distinct(trait)

rnaseq_sig_trait_either <- rnaseq_her %>%
  filter(trait %in% rnaseq_genes$trait) %>%
  filter(p_value < 0.05) %>%
  distinct(trait)

proteins_sig_trait_either <- proteins_her %>%
  filter(trait %in% proteins_genes$trait) %>%
  filter(p_value < 0.05) %>%
  distinct(trait)

genoa_sig_trait_either <- genoa_her %>%
  filter(trait %in% genoa_genes$trait) %>%
  filter(p_value < 0.05) %>%
  distinct(trait)


(nrow(rnaseq_sig_trait_either) + nrow(proteins_sig_trait_either) + nrow(genoa_sig_trait_either)) /
  (nrow(rnaseq_genes) + nrow(proteins_genes) + nrow(genoa_genes))

df_rnaseq_her <- rnaseq_her %>%
  filter(trait %in% rnaseq_genes$trait) %>%
  select(ancestry, h2g, trait) %>%
  mutate(study = "TOPMed-MESA mRNA",
    ancestry = paste0("ancestry_", ancestry)) %>%
  pivot_wider(names_from = ancestry, values_from = h2g)

df_proteins_her <- proteins_her %>%
  filter(trait %in% proteins_genes$trait) %>%
  select(ancestry, h2g, trait) %>%
  mutate(study = "TOPMed-MESA Protein",
    ancestry = paste0("ancestry_", ancestry)) %>%
  pivot_wider(names_from = ancestry, values_from = h2g)

df_genoa_her <- genoa_her %>%
  filter(trait %in% genoa_genes$trait) %>%
  select(ancestry, h2g, trait) %>%
  mutate(study = "GENOA mRNA",
    ancestry = paste0("ancestry_", ancestry)) %>%
  pivot_wider(names_from = ancestry, values_from = h2g)

p1 <- ggplot(filter(rnaseq_her, ancestry == 1 & trait %in% rnaseq_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[1]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_rnaseq_her$ancestry_1)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("EUR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p2 <- ggplot(filter(rnaseq_her, ancestry == 2  & trait %in% rnaseq_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[1]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_rnaseq_her$ancestry_2)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("AFR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p3 <- ggplot(filter(rnaseq_her, ancestry == 3  & trait %in% rnaseq_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[1]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_rnaseq_her$ancestry_3)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("HIS cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))


p4 <- ggplot(filter(proteins_her, ancestry == 1  & trait %in% proteins_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[2]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_proteins_her$ancestry_1)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("EUR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p5 <- ggplot(filter(proteins_her, ancestry == 2  & trait %in% proteins_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[2]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_proteins_her$ancestry_2)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("AFR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p6 <- ggplot(filter(proteins_her, ancestry == 3  & trait %in% proteins_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[2]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_proteins_her$ancestry_3)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("HIS cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p7 <- ggplot(filter(genoa_her, ancestry == 1 & trait %in% genoa_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[3]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_genoa_her$ancestry_1)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("EUR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p8 <- ggplot(filter(genoa_her, ancestry == 2  & trait %in% genoa_genes$trait),
  aes(x = h2g)) +
  geom_histogram(bins=100, fill = main_study_color[3]) +
  annotate("text", x = Inf, y = Inf, label = sprintf("Average cis-SNP h2g: %.2f",
    mean(df_genoa_her$ancestry_2)),
    vjust = 3, hjust = 1.2, size = 4, color = "blue") +
  ylab("Count") +
  xlab("AFR cis-SNP heritability") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
  nrow=3, ncol=3, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
  legend.grob = leg, legend = "bottom")

# ggsave("./plots/s15.png", , width = p_width+1, height = p_height+5)


slo <- tidy(lm(df_rnaseq_her$ancestry_1 ~ df_rnaseq_her$ancestry_2))[2,2]$estimate

p1 <- ggplot(df_rnaseq_her, aes(x = ancestry_1,
  y = ancestry_2)) +
  geom_point(color = main_study_color[1]) +
  xlab("EUR cis-SNP Heritability") +
  ylab("AFR cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_rnaseq_her$ancestry_1 ~ df_rnaseq_her$ancestry_3))[2,2]$estimate

p2 <- ggplot(df_rnaseq_her, aes(x = ancestry_1,
  y = ancestry_3)) +
  geom_point(color = main_study_color[1]) +
  xlab("EUR cis-SNP Heritability") +
  ylab("HIS cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_rnaseq_her$ancestry_2 ~ df_rnaseq_her$ancestry_3))[2,2]$estimate

p3 <- ggplot(df_rnaseq_her, aes(x = ancestry_1,
  y = ancestry_2)) +
  geom_point(color = main_study_color[1]) +
  xlab("AFR cis-SNP Heritability") +
  ylab("HIS cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_proteins_her$ancestry_1 ~ df_proteins_her$ancestry_2))[2,2]$estimate

p4 <- ggplot(df_proteins_her, aes(x = ancestry_1,
  y = ancestry_2)) +
  geom_point(color = main_study_color[2]) +
  xlab("EUR cis-SNP Heritability") +
  ylab("AFR cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_proteins_her$ancestry_1 ~ df_proteins_her$ancestry_3))[2,2]$estimate

p5 <- ggplot(df_proteins_her, aes(x = ancestry_1,
  y = ancestry_3)) +
  geom_point(color = main_study_color[2]) +
  xlab("EUR cis-SNP Heritability") +
  ylab("HIS cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_proteins_her$ancestry_2 ~ df_proteins_her$ancestry_3))[2,2]$estimate

p6 <- ggplot(df_proteins_her, aes(x = ancestry_2,
  y = ancestry_3)) +
  geom_point(color = main_study_color[2]) +
  xlab("AFR cis-SNP Heritability") +
  ylab("HIS cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

slo <- tidy(lm(df_genoa_her$ancestry_1 ~ df_genoa_her$ancestry_2))[2,2]$estimate

p7 <- ggplot(df_proteins_her, aes(x = ancestry_1,
  y = ancestry_2)) +
  geom_point(color = main_study_color[3]) +
  xlab("EUR cis-SNP Heritability") +
  ylab("AFR cis-SNP Heritability") +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate("text", x = 0, y = Inf, label = sprintf("Slope: %.2f", slo),
    vjust = 1.5, hjust = 0,size = 4, color = "blue") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

ggarrange(p1, p2, p3, p4, p5, p6, p7,
  nrow=3, ncol=3, labels = c("A", "B", "C", "D", "E", "F", "G"),
  legend.grob = leg, legend = "bottom")

tidy(cor.test(
  c(df_rnaseq_her$ancestry_1, df_proteins_her$ancestry_1, df_genoa_her$ancestry_1),
  c(df_rnaseq_her$ancestry_2, df_proteins_her$ancestry_2, df_genoa_her$ancestry_2)))

tidy(cor.test(
  c(df_rnaseq_her$ancestry_1, df_proteins_her$ancestry_1),
  c(df_rnaseq_her$ancestry_3, df_proteins_her$ancestry_3)))

tidy(cor.test(
  c(df_rnaseq_her$ancestry_3, df_proteins_her$ancestry_3),
  c(df_rnaseq_her$ancestry_2, df_proteins_her$ancestry_2)))

# ggsave("./plots/s16.png", width = p_width+1, height = p_height+5)

# correlation
rnaseq_corr <-
  read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_corr.tsv.gz")

proteins_corr <-
  read_tsv("~/Documents/github/data/sushie_results/real/proteins_corr.tsv.gz")

genoa_corr <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_corr.tsv.gz")

df_corr <- prepare_corr(rnaseq_corr %>%
    filter(trait %in% rnaseq_genes$trait), rnaseq_her, 3) %>%
  mutate(study = "TOPMed-MESA mRNA") %>%
  bind_rows(prepare_corr(proteins_corr %>%
      filter(trait %in% proteins_genes$trait), proteins_her, 3) %>%
      mutate(study = "TOPMed-MESA Protein"),
    prepare_corr(genoa_corr %>%
        filter(trait %in% genoa_genes$trait), genoa_her, 2) %>%
      mutate(study = "GENOA mRNA")) %>%
  mutate(name = gsub("_est_corr", "", name),
    type = factor(type, levels = c(1, 2, 3),
      labels = c("All Genes", "Sig. Heri.\n in either ancestry",
        "Sig. Heri.\n in all ancestries")),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA"))) %>%
  mutate(up_bound = corr + 1.96*se,
    low_bound = corr - 1.96*se)

pp1 <- ggplot(filter(df_corr, name == "ancestry1_ancestry2"),
  aes(x = type, y = corr, color = study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_bound, ymax = up_bound),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_y_continuous(limits = c(0.55, 1)) +
  scale_color_manual(values = main_study_color) +
  ylab("Estimated Correlations") +
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
    axis.title=element_text(face="bold", size = 8),
    axis.title.x=element_blank(),
    axis.text = element_text(face = "bold", size = 7))

pp2 <- ggplot(filter(df_corr, name == "ancestry1_ancestry3"),
  aes(x = type, y = corr, color = study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_bound, ymax = up_bound),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_y_continuous(limits = c(0.55, 1)) +
  scale_color_manual(values = main_study_color) +
  ylab("Estimated Correlations") +
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
    axis.title=element_text(face="bold", size = 8),
    axis.title.x=element_blank(),
    axis.text = element_text(face = "bold", size = 7))

pp3 <- ggplot(filter(df_corr, name == "ancestry2_ancestry3"),
  aes(x = type, y = corr, color = study)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = low_bound, ymax = up_bound),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_y_continuous(limits = c(0.55, 1)) +
  scale_color_manual(values = main_study_color) +
  ylab("Estimated Correlations") +
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
    axis.title=element_text(face="bold", size = 8),
    axis.title.x=element_blank(),
    axis.text = element_text(face = "bold", size = 7))


ggarrange(pp1, pp2, pp3, nrow=3, align = "v", label.y = 0.95, label.x = 0.05,
  labels = c("EUR vs. AFR", "EUR vs. HIS", "AFR vs. HIS"),
  font.label = list(size = 10),
  legend = "bottom", common.legend = TRUE)

# ggsave("./plots/s17.png", width = p_width-3, height = p_height+3)

# all
bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_genes$trait),
  genoa_corr %>%
    filter(trait %in% genoa_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry2_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry3_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry2_ancestry3_est_corr))


# either sig
bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait_either$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait_either$trait) %>%
    filter(trait %in% proteins_genes$trait),
  genoa_corr %>%
    filter(trait %in% genoa_sig_trait_either$trait) %>%
    filter(trait %in% genoa_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry2_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait_either$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait_either$trait) %>%
    filter(trait %in% proteins_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry3_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait_either$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait_either$trait) %>%
    filter(trait %in% genoa_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry2_ancestry3_est_corr))

# all sig
bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait$trait) %>%
    filter(trait %in% proteins_genes$trait),
  genoa_corr %>%
    filter(trait %in% genoa_sig_trait$trait) %>%
    filter(trait %in% genoa_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry2_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait$trait) %>%
    filter(trait %in% proteins_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry1_ancestry3_est_corr))

bind_rows(rnaseq_corr %>%
    filter(trait %in% rnaseq_sig_trait$trait) %>%
    filter(trait %in% rnaseq_genes$trait),
  proteins_corr %>%
    filter(trait %in% proteins_sig_trait$trait) %>%
    filter(trait %in% proteins_genes$trait)) %>%
  filter(CSIndex == 1) %>%
  summarize(m1 = mean(ancestry2_ancestry3_est_corr))

rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  distinct(trait) %>%
  mutate(sig = ifelse(trait %in% rnaseq_sig_trait$trait, 1, 0)) %>%
  bind_rows(proteins_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait) %>%
      mutate(sig = ifelse(trait %in% proteins_sig_trait$trait, 1, 0)),
    genoa_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait) %>%
      mutate(sig = ifelse(trait %in% genoa_sig_trait$trait, 1, 0))) %>%
  summarize(n = n(),
    ss = sum(sig))

9885/21088
# scatter plots
# set.seed(123)
# tmp_weight <- rnaseq_her %>%
#   filter(trait %in% rnaseq_genes$trait) %>%
#   filter(p_value < 0.05) %>%
#   group_by(trait) %>%
#   filter(n() == 3) %>%
#   distinct(trait)
# 
# tmp_weight <- tmp_weight[sample(1:nrow(tmp_weight), 100, replace = FALSE),]
# write_tsv(tmp_weight , "~/USCHPC/trash/case/rnaseq_weights.tsv", col_names = FALSE)
# 
# tmp_weight <- proteins_her %>%
#   filter(trait %in% proteins_genes$trait) %>%
#   filter(p_value < 0.05) %>%
#   group_by(trait) %>%
#   filter(n() == 3) %>%
#   distinct(trait)
# 
# tmp_weight <- tmp_weight[sample(1:nrow(tmp_weight), 100, replace = FALSE),]
# write_tsv(tmp_weight , "~/USCHPC/trash/case/proteins_weights.tsv", col_names = FALSE)
# 
# tmp_weight <- genoa_her %>%
#   filter(trait %in% genoa_genes$trait) %>%
#   filter(p_value < 0.05) %>%
#   group_by(trait) %>%
#   filter(n() == 2) %>%
#   distinct(trait)
# 
# tmp_weight <- tmp_weight[sample(1:nrow(tmp_weight), 100, replace = FALSE),]
# write_tsv(tmp_weight , "~/USCHPC/trash/case/genoa_weights.tsv", col_names = FALSE)

rnaseq_w <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_weights.tsv.gz")
proteins_w <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_weights.tsv.gz")
genoa_w <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_weights.tsv.gz")

p1 <- ggplot(rnaseq_w, aes(x = ancestry1_sushie_weight,
  y = ancestry2_sushie_weight)) +
  geom_point(color = main_study_color[1]) +
  xlab("EUR Estimated Effect Size") +
  ylab("AFR Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p2 <- ggplot(rnaseq_w, aes(x = ancestry1_sushie_weight,
  y = ancestry3_sushie_weight)) +
  geom_point(color = main_study_color[1]) +
  xlab("EUR Estimated Effect Size") +
  ylab("HIS Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p3 <- ggplot(rnaseq_w, aes(x = ancestry2_sushie_weight,
  y = ancestry3_sushie_weight)) +
  geom_point(color = main_study_color[1]) +
  xlab("AFR Estimated Effect Size") +
  ylab("HIS Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p4 <- ggplot(proteins_w, aes(x = ancestry1_sushie_weight,
  y = ancestry2_sushie_weight)) +
  geom_point(color = main_study_color[2]) +
  xlab("EUR Estimated Effect Size") +
  ylab("AFR Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p5 <- ggplot(proteins_w, aes(x = ancestry1_sushie_weight,
  y = ancestry3_sushie_weight)) +
  geom_point(color = main_study_color[2]) +
  xlab("EUR Estimated Effect Size") +
  ylab("HIS Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p6 <- ggplot(proteins_w, aes(x = ancestry2_sushie_weight,
  y = ancestry3_sushie_weight)) +
  geom_point(color = main_study_color[2]) +
  xlab("AFR Estimated Effect Size") +
  ylab("HIS Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

p7 <- ggplot(genoa_w, aes(x = ancestry1_sushie_weight,
  y = ancestry2_sushie_weight)) +
  geom_point(color = main_study_color[3]) +
  xlab("EUR Estimated Effect Size") +
  ylab("AFR Estimated Effect Size") +
  geom_abline(slope = 1, intercept = 0, color="red") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=7),
    axis.text=element_text(size = 7, face="bold"))

ggarrange(p1, p2, p3, p4, p5, p6, p7,
  nrow=3, ncol=3, labels = c("A", "B", "C", "D", "E", "F", "G"),
  legend.grob = leg, legend = "bottom")

# ggsave("./plots/s18.png", , width = p_width+1, height = p_height+5)

# correlation heterogenity

df_dens1 <- rnaseq_corr %>%
  filter(trait %in% rnaseq_genes$trait) %>%
  filter(CSIndex <=6) %>%
  select(trait, CSIndex, contains("est_corr")) %>%
  pivot_longer(cols = c(-trait, -CSIndex)) %>%
  mutate(study = "TOPMed-MESA mRNA")

df_dens2 <- proteins_corr %>%
  filter(trait %in% proteins_genes$trait) %>%
  filter(CSIndex <=6) %>%
  select(trait, CSIndex, contains("est_corr")) %>%
  pivot_longer(cols = c(-trait, -CSIndex)) %>%
  mutate(study = "TOPMed-MESA Protein")

df_dens3 <- genoa_corr %>%
  filter(trait %in% genoa_genes$trait) %>%
  filter(CSIndex <=6) %>%
  select(trait, CSIndex, contains("est_corr")) %>%
  pivot_longer(cols = c(-trait, -CSIndex)) %>%
  mutate(study = "GENOA mRNA")

df_dens <- bind_rows(df_dens1, df_dens2, df_dens3) %>%
  mutate(name = factor(name, levels = c("ancestry1_ancestry2_est_corr",
    "ancestry1_ancestry3_est_corr", "ancestry2_ancestry3_est_corr"),
    labels = c("EUR-AFR", "EUR-HIS", "AFR-HIS")),
    study = factor(study, levels = c("TOPMed-MESA mRNA",
      "TOPMed-MESA Protein", "GENOA mRNA")),
    CSIndex = factor(CSIndex, levels = 1:6, labels = paste0("L", 1:6)))

ggplot(df_dens, aes(x = value, color = study, fill = study)) +
  geom_density(alpha=0.4) +
  scale_color_manual(values = main_study_color) +
  scale_fill_manual(values = main_study_color) +
  xlab("Estimated correlation") +
  ylab("Density") +
  facet_grid(rows=vars(CSIndex), cols = vars(name), scales = "free_y") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title.x=element_text(face="bold"),
    legend.position = "bottom",
    title = element_text(size = 10, face="bold"),
    axis.text=element_text(size = 8, face="bold"))

# ggsave("./plots/s19.png", , width = p_width+1, height = p_height+5)

df_eds <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_eds.tsv")
df_pli <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_pli_new.tsv")
df_rvis <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_rvis.tsv")
df_shet <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_shet.tsv")

df_scores <- df_eds %>% 
  mutate(score = "EDS") %>%
  rename(value = EDS) %>%
  bind_rows(
    df_rvis %>%
      mutate(score = "RVIS") %>%
      rename(value = RVIS),
    df_shet %>%
      mutate(score = "s_het") %>%
      rename(value = s_het),
    df_pli %>%
      rename(score = name))

rnaseq_qtl <- rnaseq_cov %>%
  # filter(!is.na(snp)) %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  group_by(trait) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "TOPMed-MESA mRNA") %>%
  mutate(trait = gsub("_.*", "", trait)) %>%
  left_join(df_scores) %>%
  filter(!is.na(value))

proteins_qtl <- proteins_cov %>%
  # filter(!is.na(snp)) %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  group_by(trait) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "TOPMed-MESA Protein") %>%
  mutate(trait = gsub("_.*", "", trait)) %>%
  left_join(df_scores) %>%
  filter(!is.na(value))

genoa_qtl <- genoa_cov %>%
  group_by(trait) %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "GENOA mRNA") %>%
  mutate(trait = gsub("_.*", "", trait)) %>%
  left_join(df_scores) %>%
  filter(!is.na(value))

df_qtl <- bind_rows(rnaseq_qtl, proteins_qtl, genoa_qtl) %>%
  mutate(study = factor(study, levels = c("TOPMed-MESA mRNA",
    "TOPMed-MESA Protein", "GENOA mRNA")))

# constraint analysis preparation
# fst info
rnaseq_fst <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_fst.all_snp.tsv.gz") %>%
  mutate(pop1 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "EUR", `#POP1`),
    pop2 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "AFR", POP2),
    type = ifelse(pop1 == "EUR" & pop2 == "AFR", "ancestry1_ancestry2_est_corr",
      ifelse(pop1 == "EUR" & pop2 == "HIS", "ancestry1_ancestry3_est_corr",
        ifelse(pop1 == "AFR" & pop2 == "HIS", "ancestry2_ancestry3_est_corr",
          NA))),
    study = "mesa.mrna") %>%
  select(study, trait, type, all_fst = HUDSON_FST) %>%
  mutate(type = gsub("_est_corr", "", type))

# proteins
proteins_fst <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_fst.all_snp.tsv.gz") %>%
  mutate(pop1 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "EUR", `#POP1`),
    pop2 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "AFR", POP2),
    type = ifelse(pop1 == "EUR" & pop2 == "AFR", "ancestry1_ancestry2_est_corr",
      ifelse(pop1 == "EUR" & pop2 == "HIS", "ancestry1_ancestry3_est_corr",
        ifelse(pop1 == "AFR" & pop2 == "HIS", "ancestry2_ancestry3_est_corr",
          NA))),
    study = "mesa.proteins") %>%
  select(study, trait, type, all_fst = HUDSON_FST) %>%
  mutate(type = gsub("_est_corr", "", type))

# genoa
genoa_fst <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_fst.all_snp.tsv.gz") %>%
  mutate(pop1 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "EUR", `#POP1`),
    pop2 = ifelse(`#POP1` == "AFR" & POP2 == "EUR", "AFR", POP2),
    type = "ancestry1_ancestry2_est_corr",
    study = "genoa.mrna") %>%
  select(study, trait, type, all_fst = HUDSON_FST) %>%
  mutate(type = gsub("_est_corr", "", type))

all_fst <- bind_rows(rnaseq_fst, proteins_fst, genoa_fst) 


# corr info
df_rnaseq_corr <- rnaseq_corr %>%
  inner_join(rnaseq_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("corr")) %>%
  pivot_longer(cols = 3:5) %>%
  rename(type = name, corr = value) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, type, CSIndex, corr)

df_rnaseq_corr %>%
  filter(type == "ancestry1_ancestry2_est_corr") %>%
  filter(CSIndex == 1) %>%
  mutate(gene = trait,
    trait = gsub("_.+", "", trait)) %>%
  inner_join(df_scores %>%
      filter(score == "pLI") %>%
      mutate(group = ifelse(value < 0.1, "low",
        ifelse(value >0.9, "high", "middle")))) %>%
  group_by(group) %>%
  summarize(mval = mean(corr),
    sdval = sd(corr)/sqrt(n()))


df_rnaseq_corr <- rnaseq_corr %>%
  inner_join(rnaseq_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("corr")) %>%
  pivot_longer(cols = 3:5) %>%
  rename(type = name, corr = value) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, type, CSIndex, corr)


df_rnaseq_covar <- rnaseq_corr %>%
  inner_join(rnaseq_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("covar")) %>%
  pivot_longer(cols = 3:5) %>%
  rename(type = name, covar = value) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, type, CSIndex, covar)

df_proteins_corr <- proteins_corr %>%
  inner_join(proteins_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("corr")) %>%
  pivot_longer(cols = 3:5) %>%
  rename(type = name, corr = value) %>%
  mutate(study = "mesa.proteins") %>%
  select(study, trait, type, CSIndex, corr)

df_proteins_covar <- proteins_corr %>%
  inner_join(proteins_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("covar")) %>%
  pivot_longer(cols = 3:5) %>%
  rename(type = name, covar = value) %>%
  mutate(study = "mesa.proteins") %>%
  select(study, trait, type, CSIndex, covar)

df_genoa_corr <- genoa_corr %>%
  inner_join(genoa_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("corr")) %>%
  pivot_longer(cols = 3) %>%
  rename(type = name, corr = value) %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, type, CSIndex, corr)

df_genoa_covar <- genoa_corr %>%
  inner_join(genoa_cov %>%
      filter(!is.na(snp)) %>%
      distinct(trait, CSIndex)) %>%
  select(trait, CSIndex, contains("covar")) %>%
  pivot_longer(cols = 3) %>%
  rename(type = name, covar = value) %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, type, CSIndex, covar)

df_corr <- bind_rows(df_rnaseq_corr, df_proteins_corr, df_genoa_corr) %>%
  mutate(type = gsub("_est_corr", "", type))
df_covar <- bind_rows(df_rnaseq_covar, df_proteins_covar, df_genoa_covar) %>%
  mutate(type = gsub("_est_covar", "", type))

df_comp <- df_scores %>%
  inner_join(df_corr %>%
      full_join(df_covar) %>%
      mutate(gene = trait,
        trait = gsub("_.+", "", trait)), by = "trait") %>%
  left_join(all_fst %>%
      rename(gene = trait)) %>%
  mutate(CSIndex = factor(CSIndex, levels = 1:10, ordered = TRUE))

normal_comp <- df_comp %>%
  nest_by(score) %>%
  mutate(mod = list(lm(value ~ corr+CSIndex+type +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "corr")  %>%
  mutate(p.value = p.value/2) %>%
  select(score, estimate, p.value)


l1_comp <- df_comp %>%
  filter(CSIndex == 1) %>%
  nest_by(score) %>%
  mutate(mod = list(lm(value ~ corr+type +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "corr")  %>%
  mutate(p.value = p.value/2)  %>%
  select(score, estimate, p.value)

fst_comp <- df_comp %>%
  nest_by(score) %>%
  mutate(mod = list(lm(value ~ corr+CSIndex+type +study+all_fst, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "corr")  %>%
  mutate(p.value = p.value/2)  %>%
  select(score, estimate, p.value)


covar_comp <- df_comp %>%
  nest_by(score) %>%
  mutate(mod = list(lm(value ~ covar+CSIndex+type +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "covar")  %>%
  mutate(p.value = p.value/2)  %>%
  select(score, estimate, p.value)

# bootstrap
btsp_func <- function(df) {
  df_base <- tibble()
  for (cons_score in c("EDS", "RVIS", "pLI", "LOEUF", "s_het")) {
    for (idx in 1:10) {
      for (study_name in c("mesa.mrna", "mesa.proteins", "genoa.mrna")) {
        df_tmp <- df %>%
          filter(score %in% cons_score &
              CSIndex == idx & study ==  study_name)
        sel_idx <- sample(nrow(df_tmp), nrow(df_tmp), replace = TRUE)
        df_base <- df_base %>%
          bind_rows(df_tmp[sel_idx,])
      }
    }
  }
  res <- df_base %>%
    nest_by(score) %>%
    mutate(mod = list(lm(value ~ corr+CSIndex+type +study, data = data))) %>%
    reframe(tidy(mod)) %>%
    filter(term == "corr")  %>%
    select(score, estimate)
  return(res)
  
}

set.seed(123)
btsp_res <- tibble()
for (n_time in 1:100) {
  print(n_time)
  btsp_res <- btsp_res %>%
    bind_rows(btsp_func(df_comp) %>%
        mutate(times = n_time))
}

# write_tsv(btsp_res, "~/Downloads/btsp_res.tsv")

btsp_res <- read_tsv("~/Documents/github/data/sushie_results/real/btsp_res.tsv")
btsp_comp <- normal_comp %>%
  left_join(
    btsp_res %>%
      group_by(score) %>%
      summarize(sdval = sd(estimate)))  %>%
  mutate(statistic = estimate/sdval,
    p.value = pnorm(abs(statistic), lower.tail = FALSE)) %>%
  select(score, estimate, p.value)

total_comp <- normal_comp %>%
  mutate(type="Base Model") %>%
  bind_rows(btsp_comp %>%
      mutate(type="Bootstrap SE"),
    l1_comp %>%
      mutate(type="Primary Effect"),
    covar_comp %>%
      mutate(type="Effect Covariance"),
    fst_comp %>%
      mutate(type="Adjusted Fst")) %>%
  mutate(p.value = formatC(p.value, format = "e", digits = 2),
    estimate = formatC(estimate, format = "f", digits = 3),
    res = paste0(estimate," (", p.value, ")")) %>%
  select(score, type, res) %>%
  pivot_wider(names_from = score, values_from = res) %>%
  mutate(type = factor(type, levels = c("Base Model",
    "Bootstrap SE", "Primary Effect", "Effect Covariance", "Adjusted Fst"))) %>%
  select(type, pLI, LOEUF, s_het, RVIS, EDS)

write_tsv(total_comp, "~/Downloads/table1.tsv")
library(gt)

total_comp |>
  gt(rowname_col = "type") |>
  cols_label(
    s_het =  html("<b>s<sub>het</sub></b>"),
    RVIS =  html("<b>RVIS</b>"),
    EDS =  html("<b>EDS</b>"),
    pLI =  html("<b>pLI</b>"),
    LOEUF =  html("<b>LOEUF</b>"),
  ) |>
  cols_align(
    align ="center",
    columns = everything()
  ) |>
  cols_width(type ~ px(150),
    everything() ~ px(200)) |>
  gtsave("tab_1.rtf")

# heritability vs constraint

df_her_cons <- rnaseq_her %>%
  select(trait, ancestry,h2g) %>%
  mutate(study = "mesa.mrna",
    gene = trait,
    trait = gsub("_.+", "", trait)) %>%
  bind_rows(proteins_her %>%
      select(trait, ancestry,h2g) %>%
      mutate(study = "mesa.proteins",
        gene = trait,
        trait = gsub("_.+", "", trait)),
    genoa_her %>%
      select(trait, ancestry,h2g) %>%
      mutate(study = "genoa.mrna",
        gene = trait,
        trait = gsub("_.+", "", trait))) %>%
  inner_join(df_scores)

df_her_res <- df_her_cons %>%
  nest_by(score) %>%
  mutate(mod = list(lm(value ~ h2g + ancestry+ study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "h2g") %>%
  mutate(p.value = p.value /2)

# write_tsv(df_her_res, "./tables/s8.tsv")

# qtl info
rnaseq_qtl <- rnaseq_cov %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  group_by(trait) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "mesa.mrna")

proteins_qtl <- proteins_cov %>%
  # filter(!is.na(snp)) %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  group_by(trait) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "mesa.proteins") 

genoa_qtl <- genoa_cov %>%
  group_by(trait) %>%
  mutate(zero = ifelse(is.na(snp), 0, 1)) %>%
  summarize(n = ifelse(0 %in% zero, 0, length(unique(CSIndex)))) %>%
  mutate(study = "genoa.mrna")

all_qtl <- bind_rows(rnaseq_qtl, proteins_qtl, genoa_qtl)


# distance info
rnaseq_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_rnaseq_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X12, TSS = X6)

proteins_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_proteins_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X15, TSS = X6)

genoa_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/genoa_sushie_gene_list_noMHC.tsv", col_names = FALSE) %>%
  select(trait = X2, TSS = X6)

rnaseq_dist <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, pip_all, trait) %>%
  left_join(rnaseq_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  mutate(study = "mesa.mrna")

proteins_dist <- proteins_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, pip_all, trait) %>%
  left_join(proteins_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  mutate(study = "mesa.proteins")

genoa_dist <- genoa_cov %>%
  filter(!is.na(snp)) %>%
  select(snp, pos, pip_all, trait) %>%
  left_join(genoa_ref, by = "trait") %>%
  mutate(dist = abs(pos - TSS) + 1) %>%
  group_by(trait) %>%
  summarize(dist = sum(pip_all*dist) / sum(pip_all)) %>%
  mutate(study = "genoa.mrna")

all_dist <- bind_rows(rnaseq_dist, proteins_dist, genoa_dist)

# enrich score
rnaseq_enrich <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_enrich_all.tsv.gz") %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait)%>%
  filter(method == "pip_cov") %>%
  filter(converged == 1 & reg == 0) %>%
  filter(anno %in% c("dELS", "PLS", "pELS")) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, anno, est)

proteins_enrich <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_enrich_all.tsv.gz") %>%
  filter(method == "pip_cov") %>%
  filter(trait %in% filter(proteins_cov, !is.na(snp))$trait)%>%
  filter(converged == 1 & reg == 0) %>%
  filter(anno %in% c("dELS", "PLS", "pELS")) %>%
  mutate(study = "mesa.proteins") %>%
  select(study, trait, anno, est)

genoa_enrich <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_enrich_all.tsv.gz") %>%
  filter(trait %in% filter(genoa_cov, !is.na(snp))$trait)%>%
  filter(method == "pip_cov") %>%
  filter(converged == 1 & reg == 0) %>%
  filter(anno %in% c("dELS", "PLS", "pELS")) %>%
  mutate(study = "genoa.mrna") %>%
  select(study, trait, anno, est)

all_enrich <- bind_rows(rnaseq_enrich, proteins_enrich, genoa_enrich)

df_all <- all_qtl %>%
  mutate(type = "Number of molQTLs") %>%
  rename(stats = n) %>%
  bind_rows(all_dist %>%
      mutate(type = "PIP-weighted Distance (kb) to TSS",
        dist = dist/1000) %>%
      rename(stats = dist),
    all_enrich %>%
      mutate(anno = paste0(anno, " Log Enrichment")) %>%
      rename(stats = est,
        type = anno)) %>%
  mutate(gene = trait,
    trait = gsub("_.+", "", trait)) %>%
  inner_join(df_scores) %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF",
    "s_het", "RVIS", "EDS")),
    type = factor(type, levels = c("Number of molQTLs",
      "PIP-weighted Distance (kb) to TSS", "PLS Log Enrichment",
      "pELS Log Enrichment", "dELS Log Enrichment"))) %>%
  left_join(all_fst %>%
      rename(gene = trait) %>%
      distinct(study, gene, all_fst))

btsp_func2 <- function(df) {
  df_base <- tibble()
  for (cons_score in c("EDS", "RVIS", "pLI", "LOEUF", "s_het")) {
    for (want_type in c("Number of molQTLs",
      "PIP-weighted Distance (kb) to TSS", "pELS Log Enrichment",
      "PLS Log Enrichment", "dELS Log Enrichment")) {
      for (study_name in c("mesa.mrna", "mesa.proteins", "genoa.mrna")) {
        df_tmp <- df %>%
          filter(score %in% cons_score &
              type == want_type & study ==  study_name)
        sel_idx <- sample(nrow(df_tmp), nrow(df_tmp), replace = TRUE)
        df_base <- df_base %>%
          bind_rows(df_tmp[sel_idx,])
      }
    }
  }
  res <- df_base %>%
    nest_by(type, score) %>%
    mutate(mod = list(lm(value ~ stats +study, data = data))) %>%
    reframe(tidy(mod)) %>%
    filter(term == "stats")  %>%
    select(type, score, estimate)
  return(res)
}

set.seed(123)
btsp_res2 <- tibble()
for (n_time in 1:100) {
  print(n_time)
  btsp_res2 <- btsp_res2 %>%
    bind_rows(btsp_func2(df_all) %>%
        mutate(times = n_time))
}
# write_tsv(btsp_res2, "~/Downloads/btsp_res2.tsv")

btsp_res2 <- read_tsv("~/Documents/github/data/sushie_results/real/btsp_res2.tsv")

df_text <- df_all %>%
  nest_by(type, score) %>%
  mutate(mod = list(lm(value ~ stats +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "stats")  %>%
  left_join(btsp_res2 %>%
      group_by(type, score) %>%
      summarize(sdval = sd(estimate))) %>%
  left_join(df_all %>%
      nest_by(type, score) %>%
      mutate(mod = list(lm(value ~ stats +study+all_fst, data = data))) %>%
      reframe(tidy(mod)) %>%
      filter(term == "stats") %>%
      select(type, score, fst_p = p.value)) %>%
  mutate(p.value = pnorm(abs(estimate/sdval), lower.tail = FALSE),
    label = sprintf("beta = %.2e\np = %.2e\nFst-adjusted p=%.2e", estimate, p.value,fst_p)) %>%
  select(type, score, label)

df_n_qtl <- df_all %>%
  filter(type == "Number of molQTLs") %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF",
    "s_het", "RVIS", "EDS")))

df_n_qtl_text <- df_text %>%
  filter(type == "Number of molQTLs") %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF",
    "s_het", "RVIS", "EDS")))

ggplot(df_n_qtl, aes(x = stats, y = value)) +
  geom_bin_2d(bins=200) +
  geom_smooth(method = "lm") +
  facet_wrap(~score, ncol=2, scales = "free") +
  ylab("Constraint Score") +
  xlab("Number of molQTLs") +
  scale_x_continuous(breaks=0:10) +
  geom_text(data = df_n_qtl_text, aes(label = label, x = Inf, y = Inf), 
    vjust = 1, hjust = 1, size = 1.5, color = "red") +
  scale_fill_continuous(name = "Density") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.text = element_text(size = 5, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 4),
    title = element_text(size = 5, face="bold"),
    legend.key.size = unit(0.4, 'cm'),
    axis.text=element_text(size = 5, face="bold")) 

# ggsave("./plots/s20.png", width = p_width/1.75, height = p_height+1)

df_n_qtl <- df_all %>%
  filter(type == "PIP-weighted Distance (kb) to TSS") %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF",
    "s_het", "RVIS", "EDS")),
    type = "Expected Distance (kb) to TSS")

df_n_qtl_text <- df_text %>%
  filter(type == "PIP-weighted Distance (kb) to TSS") %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF",
    "s_het", "RVIS", "EDS")),
    type = "Expected Distance (kb) to TSS")

ggplot(df_n_qtl, aes(x = stats, y = value)) +
  geom_bin_2d(bins=200) +
  geom_smooth(method = "lm") +
  facet_wrap(~score, ncol=2, scales = "free") +
  ylab("Constraint Score") +
  xlab("PIP-weighted Distance (kb) to TSS") +
  geom_text(data = df_n_qtl_text, aes(label = label, x = Inf, y = Inf), 
    vjust = 1, hjust = 1, size = 1.5, color = "red") +
  scale_fill_continuous(name = "Density") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.text = element_text(size = 5, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 4),
    title = element_text(size = 5, face="bold"),
    legend.key.size = unit(0.4, 'cm'),
    axis.text=element_text(size = 5, face="bold")) 

# ggsave("./plots/s21.png", width = p_width/1.75, height = p_height+1)


df_n_qtl <- df_all %>%
  filter(grepl("Enrich", type)) %>%
  mutate(type = factor(type, levels = c("PLS Log Enrichment",
    "pELS Log Enrichment", "dELS Log Enrichment")),
    score = factor(score, levels = c("pLI", "LOEUF",
      "s_het", "RVIS", "EDS")))

df_n_qtl_text <- df_text %>%
  filter(grepl("Enrich", type)) %>%
  mutate(type = factor(type, levels = c("PLS Log Enrichment",
    "pELS Log Enrichment", "dELS Log Enrichment")),
    score = factor(score, levels = c("pLI", "LOEUF",
      "s_het", "RVIS", "EDS")))

ggplot(df_n_qtl, aes(x = stats, y = value)) +
  geom_bin_2d(bins=200) +
  geom_smooth(method = "lm") +
  facet_grid(rows=vars(score),cols=vars(type), scales = "free") +
  ylab("Constraint Score") +
  geom_text(data = df_n_qtl_text, aes(label = label, x = Inf, y = Inf), 
    vjust = 1, hjust = 1, size = 1.5, color = "red") +
  scale_fill_continuous(name = "Density") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.position = "bottom",
    strip.text = element_text(size = 5, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    axis.title.x=element_blank(),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 4),
    title = element_text(size = 5, face="bold"),
    legend.key.size = unit(0.4, 'cm'),
    axis.text=element_text(size = 5, face="bold")) 

# ggsave("./plots/s22.png", width = p_width/1.75, height = p_height+2)

