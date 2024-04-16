library(tidyverse)
library(ggpubr)
library(broom)
source("./utils.R")

rnaseq_cov <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz")

proteins_cov <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz")

genoa_cov <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz")

twas_colors <- c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#a6cee3", "SuSiE" = "#e7298a",
  "LASSO" = "#66a61e", "Elastic Net" = "#e6ab02", "gBLUP" = "#a6761d")

rnaseq_r2 <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_r2.tsv.gz") %>%
  mutate(study = "mesa.mrna") %>%
  filter(trait %in% filter(rnaseq_cov, !is.na(snp))$trait)

proteins_r2 <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_r2.tsv.gz") %>%
  mutate(study = "mesa.proteins") %>%
  filter(trait %in% filter(proteins_cov, !is.na(snp))$trait)

genoa_r2 <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_r2.tsv.gz") %>%
  mutate(study = "genoa.mrna") %>%
  filter(trait %in% filter(genoa_cov, !is.na(snp))$trait)

rnaseq_corr <-
  read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_corr.tsv.gz")

proteins_corr <-
  read_tsv("~/Documents/github/data/sushie_results/real/proteins_corr.tsv.gz")

genoa_corr <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_corr.tsv.gz")


df_r2 <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(type == "r2") %>%
  select(-cross) %>%
  pivot_longer(cols = 1:7) %>%
  group_by(name) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    study = factor(study,
      levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Proteins",
        "GENOA mRNA")))

r2_res <- tibble()
for (met in c("SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")) {
  df_tmp <- df_r2 %>%
    filter(name %in% c("SuShiE", met)) %>%
    mutate(name = factor(name, levels = c("SuShiE", met)))
  
  r2_res <- r2_res %>%
    bind_rows(tidy(lm(value ~ name + study, df_tmp))[2,])
}

df1 <- r2_res %>%
  mutate(p.value = ifelse(estimate < 0, pnorm(abs(statistic), lower.tail = FALSE),
    pnorm(-abs(statistic), lower.tail = FALSE)))


heter_genes <- bind_rows(
  rnaseq_corr %>%
    inner_join(rnaseq_cov %>%
        filter(!is.na(snp)) %>%
        distinct(trait, CSIndex)) %>%
    select(trait, CSIndex, contains("corr")) %>%
    pivot_longer(cols = 3:5) %>%
    rename(type = name, corr = value) %>%
    mutate(study = "mesa.mrna") %>%
    select(study, trait, type, CSIndex, corr) %>%
    filter(corr < 0.5) %>%
    distinct(trait),
  proteins_corr %>%
    inner_join(proteins_cov %>%
        filter(!is.na(snp)) %>%
        distinct(trait, CSIndex)) %>%
    select(trait, CSIndex, contains("corr")) %>%
    pivot_longer(cols = 3:5) %>%
    rename(type = name, corr = value) %>%
    mutate(study = "mesa.mrna") %>%
    select(study, trait, type, CSIndex, corr) %>%
    filter(corr < 0.5) %>%
    distinct(trait),
  genoa_corr %>%
    inner_join(genoa_cov %>%
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

df_r2_heter <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(type == "r2") %>%
  select(-cross) %>%
  pivot_longer(cols = 1:7) %>%
  group_by(name) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    study = factor(study,
      levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Proteins",
        "GENOA mRNA"))) %>%
  filter(trait %in% heter_genes$trait)

r2_res_heter <- tibble()
for (met in c("SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")) {
  df_tmp <- df_r2_heter %>%
    filter(name %in% c("SuShiE", met)) %>%
    mutate(name = factor(name, levels = c("SuShiE", met)))
  
  r2_res_heter <- r2_res_heter %>%
    bind_rows(tidy(lm(value ~ name + study, df_tmp))[2,])
}

df2 <- r2_res_heter %>%
  mutate(p.value = ifelse(estimate < 0, pnorm(abs(statistic), lower.tail = FALSE),
    pnorm(-abs(statistic), lower.tail = FALSE)))

df_r2_comp <- df_r2 %>%
  group_by(name) %>%
  summarize(mval = mean(value),
    se_upp = mval + 1.96 * sd(value)/sqrt(n()),
    se_low = mval - 1.96 * sd(value)/sqrt(n())) %>%
  mutate(type = "all") %>%
  bind_rows(df_r2_heter %>%
      group_by(name) %>%
      summarize(mval = mean(value),
        se_upp = mval + 1.96 * sd(value)/sqrt(n()),
        se_low = mval - 1.96 * sd(value)/sqrt(n())) %>%
      mutate(type = "heter")) %>%
  mutate(type = factor(type, levels = c("all", "heter"),
    labels = c("A: All e/pGenes", "B:e/pGenes Exhibited Heterogeneity")))

ggplot(df_r2_comp,
  aes(x = name, y = mval, color = name)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = se_low, ymax = se_upp),
    position=position_dodge(width=0.5), width = 0.2) +
  facet_grid(cols = vars(type)) +
  scale_color_manual(values = twas_colors) +
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
    title = element_text(size = 10, face="bold"),
    axis.text=element_text(size = 8, face="bold"))

# ggsave("./plots/s23.png", width = p_width, height = p_height)

df_r2_cross <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(type == "r2") %>%
  select(sushie, cross, type, trait, study) %>%
  pivot_longer(cols = 1:2) %>%
  group_by(name) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(name = factor(name,
    levels = c("sushie", "cross"),
    labels = c("Ancestry-matched weights", "Cross-ancestry weights")),
    study = factor(study,
      levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Proteins",
        "GENOA mRNA"))) %>%
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

# ggsave("./plots/s24.png", width = p_width/2, height = p_height)

df_r2_cross2 <- bind_rows(rnaseq_r2,
  proteins_r2,
  genoa_r2) %>%
  filter(type == "r2") %>%
  select(sushie, cross, type, trait, study) %>%
  pivot_longer(cols = 1:2) %>%
  group_by(name) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(name = factor(name,
    levels = c("sushie", "cross"),
    labels = c("Ancestry-matched weights", "Cross-ancestry weights")),
    study = factor(study,
      levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Proteins",
        "GENOA mRNA")))

tidy(lm(value ~ name + study, df_r2_cross2))

rnaseq_cov_simple <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.mrna") %>%
  distinct()

proteins_cov_simple <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz")%>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.proteins") %>%
  distinct()

genoa_cov_simple <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "genoa.mrna") %>%
  distinct()

cov_genes_list <- bind_rows(rnaseq_cov_simple,
  proteins_cov_simple,
  genoa_cov_simple)

nrow(rnaseq_cov_simple) + nrow(genoa_cov_simple)
nrow(proteins_cov_simple)


rnaseq_mega_simple <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.mega_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.mrna") %>%
  distinct()

proteins_mega_simple <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.mega_cs.tsv.gz")%>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.proteins") %>%
  distinct()

genoa_mega_simple <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.mega_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "genoa.mrna") %>%
  distinct()

mega_genes_list <- bind_rows(rnaseq_mega_simple,
  proteins_mega_simple,
  genoa_mega_simple)


wbc_traits <- c("WBC", "MON", "NEU", "EOS", "BAS", "LYM")
df_twas <- read_tsv(
  "~/Documents/github/data/sushie_results/real/sushie_twas.tsv.gz") %>%
  filter(pheno %in% wbc_traits)

df_twas %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur + n_afr + n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  summarize(m = mean(n_total))

df_twas %>%
  inner_join(cov_genes_list) %>%
  pivot_longer(cols = c(n_eur, n_afr, n_his)) %>%
  select(pheno, study, name, value) %>%
  filter(study != "genoa.mrna") %>%
  group_by(name) %>%
  summarize(m = mean(value))

df_twas %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur+n_afr+n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  group_by(pheno) %>%
  summarize(m = mean(n_total))

df_twas %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur+n_afr+n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  # group_by(pheno) %>%
  summarize(m = mean(n_total))

104242 - mean(c(79589, 82238, 84805, 80885, 86313))

cov_twas <- df_twas %>%
  filter(term == "df_anc") %>%
  inner_join(cov_genes_list) 

cov_twas %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

cov_twas %>%
  filter(p.value < 0.05/23000) %>%
  group_by(pheno) %>%
  summarize(n = n())

199/221

sig_res <- cov_twas %>%
  filter(p.value < 0.05/23000)  %>%
  select(pheno, study, gene, estimate, statistic, p.value, n_eur, n_afr, n_his)

# write_tsv(sig_res, "./tables/s10.tsv")


# lu et al
df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/gencode.v34.gene.only.tsv.gz")

load("~/Documents/github/MA-FOCUS-data-code/real-data/data/twas.RData")

lu_twas_sig <- twas_all %>%
  filter(PHEN %in% wbc_traits) %>%
  filter(TWAS.P < 0.05/23000) %>%
  select(gene = ID, pheno = PHEN, twas = TWAS.Z, POP) %>%
  left_join(df_ref %>%
      select(gene = NAME, ID2)) %>%
  filter(!is.na(ID2)) %>%
  select(-gene) %>%
  rename(gene = ID2) %>%
  distinct(pheno, gene)

lu_twas_sig %>%
  group_by(pheno) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

sig_thresh <- qnorm(0.05/23000, lower.tail = FALSE)

# kachuri et al
lfs <- list.files("~/Downloads/Results/", pattern="Result", full.names = TRUE)
kachuri_twas_sig <- lfs %>% map_df(read_csv) %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "WBC") %>%
  filter(abs(zscore) > sig_thresh) %>%
  distinct(gene, pheno)

mon_twas <- read_csv("~/Downloads/Results/MONO_AA_whole_blood.csv") %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "MON") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

neu_twas <- read_csv("~/Downloads/Results/NEU_AA_whole_blood.csv") %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "NEU") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

wbc_twas <- read_csv("~/Downloads/Results/WBC_AA_whole_blood.csv") %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "WBC") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

kachuri_twas_sig2 <- bind_rows(mon_twas, neu_twas, wbc_twas)

# tapia et al 
library(readxl)
tapia_twas <- read_excel("~/Downloads/gepi22436-sup-0003-supptable2_twassignif.xlsx",sheet = 1, skip = 3)

# lymphocyte count
# monocyte count
# neutrophil count
# white blood cell count

tapia_twas_sig <- tapia_twas %>%
  mutate(pheno = ifelse(Trait %in% "lymphocyte count", "LYM",
    ifelse(Trait %in% "monocyte count", "MON",
      ifelse(Trait %in% "neutrophil count", "NEU",
        ifelse(Trait %in% "white blood cell count", "WBC", NA))))) %>%
  filter(!is.na(pheno)) %>%
  filter(`p...20`< 2.173913e-06) %>%
  distinct(Gene, pheno) %>%
  select(NAME = Gene, pheno) %>%
  left_join(df_ref %>%
      select(NAME, gene = ID2)) %>%
  select(gene, pheno)

tapia_twas_sig %>%
  group_by(pheno) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# rowland et al 

rowland_twas <- read_excel("~/Downloads/supplemental_tables_ddac011.xlsx",sheet = 2) %>%
  filter(marginal_significant == 1) %>%
  filter(pheno_class == "WBC")

# "basophil_count"        
# "eosinophil_count"      
# "lymphocyte_count"      
# "monocyte_count"        
# "neutrophil_count"      
# "white_blood_cell_count"

rowland_twas_sig <- rowland_twas %>%
  mutate(pheno = ifelse(phenotype %in% "basophil_count", "BAS",
    ifelse(phenotype %in% "eosinophil_count", "EOS",
      ifelse(phenotype %in% "lymphocyte_count", "LYM",
        ifelse(phenotype %in% "monocyte_count", "MON",
          ifelse(phenotype %in% "neutrophil_count", "NEU",
            ifelse(phenotype %in% "white_blood_cell_count", "WBC", NA))))))) %>%
  left_join(df_ref %>%
      select(gene = NAME, ID2), by = c("gene_name"="gene")) %>%
  select(gene= ID2, pheno) %>%
  filter(!is.na(gene))

rowland_twas_sig %>%
  group_by(pheno) %>%
  summarize(n = n()) %>%
  arrange(desc(n))

# wen et al
library(stringr)
wen_twas_sig <- read_excel("~/Downloads/genes-12-01049-s001 (1)/Supplementary_Tables_finalSUB_R2.xlsx", sheet = 4, skip = 1) %>%
  mutate(gene = str_extract(ENSGID, "ENSG.+")) %>%
  select(gene, pheno = Phenotype) %>%
  filter(pheno %in% "WBC")

all_twas_sig <- bind_rows(
  lu_twas_sig,
  kachuri_twas_sig,
  kachuri_twas_sig2,
  tapia_twas_sig,
  rowland_twas_sig,
  wen_twas_sig
) %>%
  distinct() 

all_twas_sig %>%
  group_by(pheno) %>%
  summarize(n = n())

mega_twas <- df_twas %>%
  filter(term == "df_meg") %>%
  inner_join(mega_genes_list)

mean(cov_twas$statistic^2)
mean(mega_twas$statistic^2)


mega_twas %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

221-177
221/177

cov_twas %>%
  filter(study != "mesa.proteins") %>%
  filter(p.value < 0.05/23000) %>%
  select(gene, study, pheno) %>%
  mutate(full_gene = gene,
    gene = gsub("_.+", "", gene)) %>%
  left_join(all_twas_sig %>%
      mutate(rep = 1)) %>%
  summarize(num = sum(!is.na(rep)),
    dem = n(),
    ratio = num/dem)

mega_twas %>%
  filter(study != "mesa.proteins") %>%
  filter(p.value < 0.05/23000) %>%
  select(gene, study, pheno) %>%
  mutate(full_gene = gene,
    gene = gsub("_.+", "", gene)) %>%
  left_join(all_twas_sig %>%
      mutate(rep = 1)) %>%
  summarize(num = sum(!is.na(rep)),
    dem = n(),
    ratio = num/dem)

prop.test(c(119, 101), c(211, 169))

sig_thresh <- qnorm(0.05/23000, lower.tail = FALSE)

df_twas <- cov_twas %>%
  select(study, gene, pheno, sushie = statistic) %>%
  inner_join(mega_twas %>%
      select(study, gene, pheno, susie = statistic)) %>%
  mutate(type = ifelse(abs(sushie) > sig_thresh &
      abs(susie) > sig_thresh, "Both",
    ifelse(abs(sushie) > sig_thresh, "SuShiE",
      ifelse(abs(susie) > sig_thresh, "SuSiE",
        "Neither"))),
    type = factor(type, levels = c("Both",
      "SuShiE",
      "SuSiE",
      "Neither")))

p1 <- ggplot(df_twas, aes(x = susie, y = sushie, color = type)) +
  geom_point() +
  geom_smooth(method = "lm", color="#f4a582", alpha=0.5) +
  scale_color_manual(values = c("#a6611a", "#1b9e77", "#e7298a", "lightgrey")) +
  xlab(expression(bold("SuSiE T/PWAS" ~ t))) +
  ylab(expression(bold("SuShiE T/PWAS" ~ t))) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 8),
    legend.position = "bottom",
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 7),
    axis.text = element_text(face = "bold", size = 6)) 


df_eds <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_eds.tsv")
df_pli <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_pli_new.tsv")
df_rvis <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_rvis.tsv")
df_shet <- read_tsv("~/Documents/github/data/sushie_results/Constraint/df_shet.tsv")

df_scores <- df_eds %>% 
  mutate(score = "EDS",
    q25 = quantile(EDS, 0.1),
    q75 = quantile(EDS, 0.9),
    Cate = ifelse(EDS < q25, "Low",
      ifelse(EDS > q75, "High", "Middle"))) %>%
  rename(value = EDS) %>%
  select(trait, value, score, Cate) %>%
  bind_rows(
    df_rvis %>% 
      mutate(score = "RVIS",
        q25 = quantile(RVIS, 0.1),
        q75 = quantile(RVIS, 0.9),
        Cate = ifelse(RVIS < q25, "High",
          ifelse(RVIS > q75, "Low", "Middle"))) %>%
      rename(value = RVIS) %>%
      select(trait, value, score, Cate) ,
    df_shet %>% 
      mutate(score = "s_het",
        q25 = quantile(s_het, 0.1),
        q75 = quantile(s_het, 0.9),
        Cate = ifelse(s_het < q25, "Low",
          ifelse(s_het > q75, "High", "Middle"))) %>%
      rename(value = `s_het`) %>%
      select(trait, value, score, Cate), 
    df_pli %>%
      rename(score = name) %>%
      group_by(score) %>%
      mutate(q25 = quantile(value, 0.1),
        q75 = quantile(value, 0.9)) %>%
      ungroup() %>%
      mutate(Cate = ifelse(score == "pLI",
        ifelse(value < 0.1, "Low",
          ifelse(value > 0.9, "High", "Middle")),
        ifelse(value < q25, "High",
          ifelse(value > q75, "Low", "Middle"))))%>%
      select(trait, value, score, Cate))



df_wk1 <- cov_twas %>%
  select(type = term, statistic, gene, pheno, study) %>%
  mutate(chisq_twas = statistic^2,
    trait = gene,
    gene = gsub("_.*", "", gene)) %>%
  left_join(df_scores %>%
      rename(gene = trait)) %>%
  filter(!is.na(value))

df_wk2 <- mega_twas %>%
  select(type = term, statistic, gene, pheno, study) %>%
  mutate(chisq_twas = statistic^2,
    trait = gene,
    gene = gsub("_.*", "", gene)) %>%
  left_join(df_scores %>%
      rename(gene = trait)) %>%
  filter(!is.na(value))

df_wk <- bind_rows(df_wk1, df_wk2)


df_asso <- df_wk %>%
  nest_by(type, score) %>%
  mutate(mod = list(lm(value ~ chisq_twas +pheno +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "chisq_twas")  %>%
  mutate(p.value = p.value/2,
    score = factor(score, levels = c("pLI", "LOEUF", "s_het", "RVIS", "EDS")),
    type = factor(type, levels = c("df_anc", "df_meg"))) %>%
  arrange(type, score)

pnorm(abs((sum(df_asso$statistic[1:5]^2) - sum(df_asso$statistic[6:10]^2))/sqrt(20)),
  lower.tail = FALSE)


# write_tsv(df_asso, "./tables/s11.tsv")


tt1 <- df_wk1 %>%
  group_by(Cate, score) %>%
  summarize(mvalue = mean(chisq_twas),
    se = sd(chisq_twas)/sqrt(n()),
    low_bound = mvalue - 1.96 *se,
    upp_bound = mvalue+1.96*se) %>%
  mutate(method = "SuShiE")

tt2 <- df_wk2 %>%
  group_by(Cate, score) %>%
  summarize(mvalue = mean(chisq_twas),
    se = sd(chisq_twas)/sqrt(n()),
    low_bound = mvalue - 1.96 *se,
    upp_bound = mvalue+1.96*se) %>%
  mutate(method = "SuSiE")


tt <- bind_rows(tt1, tt2) %>%
  mutate(method = factor(method, levels = c("SuShiE", "SuSiE")),
    score = factor(score, levels = c("pLI", "LOEUF", "s_het", "RVIS", "EDS"),
      labels = c("bold(pLI)", "bold(LOEUF)", "bold(s[het])", "bold(RVIS)",
        "bold(EDS)")),
    Cate = factor(Cate, levels = c("Low", "Middle", "High")))


p2 <- ggplot(tt, aes(x = Cate, y = mvalue, fill = method, color = method)) +
  geom_point(position=position_dodge(width=0.5), shape=21, size=2) +
  geom_errorbar(aes(ymin = low_bound, ymax = upp_bound),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_fill_manual(values = c("#1b9e77", "#e7298a")) +
  scale_color_manual(values = c("#1b9e77", "#e7298a")) +
  ylab(expression(bold("Average T/PWAS "~ Chi^2))) +
  xlab("Constraint Group") +
  facet_grid(cols = vars(score), labeller = label_parsed) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 6),
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    # legend.box.background = element_rect(colour = "black"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size = 7),
    axis.text = element_text(face = "bold", size = 6)) 
library(patchwork)
p1+p2 +plot_layout(design = "AABBBBB", guides = "collect") &
  theme(legend.position="bottom")
# + plot_layout(guides = 'collect')
# ggarrange(p1, p2, common.legend = TRUE, align = "h", legend = "bottom", labels = c("A", "B"), widths = c(1, 2.5))

# ggsave("./plots/p5.png", width = p_width, height = p_height+0.7)


df_shet_new <- df_shet %>%
  mutate(quant = ntile(s_het, 10)) %>%
  filter(quant %in% c(1, 10))

haha <- cov_twas %>%
  mutate(chisq = statistic^2) %>%
  select(term, chisq, gene, pheno) %>%
  bind_rows(mega_twas %>%
      mutate(chisq = statistic^2) %>%
      select(term, chisq, gene, pheno)) %>%
  mutate(trait = gsub("_.+", "", gene)) %>%
  filter(trait %in% filter(df_shet_new , quant ==10)$trait) %>%
  pivot_wider(names_from = term, values_from = chisq) %>%
  mutate(diff = df_anc - df_meg) %>%
  filter(!is.na(diff))


# qqplot
total_twas <- bind_rows(cov_twas, mega_twas) %>%
  mutate(term = factor(term, levels = c("df_anc", "df_meg"),
    labels = c("SuShiE T/PWAS", "SuSiE T/PWAS")),
    study = factor(study, levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
      labels = c("TOPMed-MESA mRNA",
        "TOPMed-MESA Protein", "GENOA mRNA")))

ggplot(total_twas, aes(sample = -log(p.value))) +
  stat_qq(aes(color = term), distribution = qexp) +
  stat_qq_line(distribution = qexp) +
  scale_color_manual(values = c("#1b9e77", "#e7298a")) +
  facet_grid(rows = vars(pheno), cols = vars(study), scales = "free") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 6),
    legend.position = "bottom",
    # legend.box.background = element_rect(colour = "black"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_blank(),
    axis.text = element_text(face = "bold", size = 6)) 

# ggsave("./plots/s21.png", width = p_width-1, height = p_height+3)



