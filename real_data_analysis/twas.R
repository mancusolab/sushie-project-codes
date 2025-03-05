library(tidyverse)
library(ggpubr)
library(broom)
source("./utils.R")

# change the data folder to the zenodo-downloaded data folder
data_folder <- "~/Downloads/sushie_real_data_results/all_results/"
metadata_folder <- "~/Downloads/sushie_real_data_results/metadata/"
constraint_folder <- "~/Downloads/sushie_real_data_results/constraint_score/"
validation_folder <- "~/Downloads/sushie_real_data_results/twas_validation/"

rnaseq_cov <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz"))

proteins_cov <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz"))

genoa_cov <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz"))

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

twas_colors <- c(method_colors,
  "LASSO" = "#fb8072", "Elastic Net" = "#b3de69", "gBLUP" = "#fccde5")

# lfs <- list.files("~/Downloads/new_res/", full.names = TRUE)
# 
# df_twas <- lfs %>% map_df(read_tsv, col_type = cols())
# write_tsv(df_twas, "{data_folder}/sushie_twas.tsv.gz")

rnaseq_cov_simple <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.mrna") %>%
  distinct()

proteins_cov_simple <- read_tsv(glue("{data_folder}/proteins_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.proteins") %>%
  distinct()

genoa_cov_simple <- read_tsv(glue("{data_folder}/genoa_normal.sushie_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "genoa.mrna") %>%
  distinct()

cov_genes_list <- bind_rows(rnaseq_cov_simple,
  proteins_cov_simple,
  genoa_cov_simple)

nrow(rnaseq_cov_simple) + nrow(genoa_cov_simple)
nrow(proteins_cov_simple)

rnaseq_mega_simple <- read_tsv(glue("{data_folder}/rnaseq_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.mrna") %>%
  distinct()

proteins_mega_simple <- read_tsv(glue("{data_folder}/proteins_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "mesa.proteins") %>%
  distinct()

genoa_mega_simple <- read_tsv(glue("{data_folder}/genoa_normal.mega_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  select(gene = trait) %>%
  mutate(study = "genoa.mrna") %>%
  distinct()

mega_genes_list <- bind_rows(rnaseq_mega_simple,
  proteins_mega_simple,
  genoa_mega_simple)

mesusie_genes_list <- bind_rows(
  read_tsv(glue("{data_folder}/rnaseq_mesusie_cs.tsv.gz")) %>%
    distinct(trait) %>%
    rename(gene = trait) %>%
    mutate(study = "mesa.mrna"),
  read_tsv(glue("{data_folder}/proteins_mesusie_cs.tsv.gz")) %>%
    distinct(trait) %>%
    rename(gene = trait) %>%
    mutate(study = "mesa.proteins"),
  read_tsv(glue("{data_folder}/genoa_mesusie_cs.tsv.gz")) %>%
    distinct(trait) %>%
    rename(gene = trait) %>%
    mutate(study = "genoa.mrna")
)

df_twas <- read_tsv(glue("{data_folder}/sushie_twas.tsv.gz"))

wbc_traits <- c("WBC", "MON", "NEU", "EOS", "BAS", "LYM")

twas_out <- df_twas %>%
  filter(term %in% "sushie") %>%
  filter(pheno %in% wbc_traits) %>%
  filter(p.value <= 0.05/23000) %>%
  select(Phenotype = pheno,
    Study = study,
    Gene = gene,
    `Estimate (Beta)` = estimate,
    Statistic = statistic,
    `Two-sided P-value` = p.value,
    `EUR Sample Size` = n_eur,
    `AFR Sample Size` = n_afr,
    `HIS Sample Size` = n_his) %>%
  mutate(Study = factor(Study, 
    levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
    labels = c("TOPMed-MESA mRNA", "TOPMed-MESA Proteins", "GENOA mRNA"))) %>%
  arrange(Phenotype, Study, Gene)

# write_tsv(twas_out, "./tables/s14.tsv")

# qqplot
qq_twas <-  df_twas %>%
  filter(term %in% "sushie") %>%
  filter(pheno %in% wbc_traits) %>%
  mutate(term = "SuShiE") %>%
  mutate(study = factor(study, 
    levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
    labels = c("TOPMed-MESA mRNA", "TOPMed-MESA Proteins", "GENOA mRNA")))

ggplot(qq_twas, aes(sample = -log(p.value))) +
  stat_qq(aes(color = term), distribution = qexp) +
  stat_qq_line(distribution = qexp) +
  scale_color_manual(values = twas_colors) +
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

# ggsave("./plots/s27.png", width = p_width-1, height = p_height+3)

df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur + n_afr + n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  summarize(m = mean(n_total))

df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  inner_join(cov_genes_list) %>%
  pivot_longer(cols = c(n_eur, n_afr, n_his)) %>%
  select(pheno, study, name, value) %>%
  filter(study != "genoa.mrna") %>%
  group_by(name) %>%
  summarize(m = mean(value))

df_twas %>%
  filter(!pheno %in% wbc_traits) %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur + n_afr + n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  group_by(pheno) %>%
  summarize(m = mean(n_total))

df_twas %>%
  filter(!pheno %in% wbc_traits) %>%
  inner_join(cov_genes_list) %>%
  pivot_longer(cols = c(n_eur, n_afr, n_his)) %>%
  select(pheno, study, name, value) %>%
  filter(study != "genoa.mrna") %>%
  group_by(name, pheno) %>%
  summarize(m = mean(value))

df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  inner_join(cov_genes_list) %>%
  mutate(n_total = n_eur+n_afr+n_his) %>%
  select(pheno, study, n_total) %>%
  filter(study != "genoa.mrna") %>%
  group_by(pheno) %>%
  summarize(m = mean(n_total))

104242 - mean(c(79589, 82238, 84805, 80885, 86313))

cov_twas <- df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  filter(term == "sushie") %>%
  inner_join(cov_genes_list) 

cov_twas %>%
  filter(pheno %in% wbc_traits) %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

cov_twas %>%
  filter(pheno %in% wbc_traits) %>%
  filter(p.value < 0.05/23000) %>%
  group_by(pheno) %>%
  summarize(n = n())

179/195

mega_twas <- df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  filter(term == "mega") %>%
  inner_join(mega_genes_list)

mesusie_twas <- df_twas %>%
  filter(pheno %in% wbc_traits) %>%
  filter(term == "mesusie") %>%
  inner_join(mesusie_genes_list)

df_wbc_tmp1 <- df_twas %>%
  filter(term %in% c("sushie", "mega")) %>%
  mutate(term = factor(term, levels = c("sushie", "mega"))) %>%
  filter(pheno %in% wbc_traits) %>%
  group_by(pheno, gene) %>%
  filter(n() == 2) %>%
  mutate(s2 = statistic^2) %>%
  select(term, s2, gene, pheno, study)

df_wbc_tmp2 <- df_twas %>%
  filter(term %in% c("sushie", "mesusie")) %>%
  mutate(term = factor(term, levels = c("sushie", "mesusie"))) %>%
  filter(pheno %in% wbc_traits) %>%
  group_by(pheno, gene) %>%
  filter(n() == 2) %>%
  mutate(s2 = statistic^2) %>%
  select(term, s2, gene, pheno, study)

# perform bootstrap for the standard error
set.seed(123)
n_rep <- 100
c1_comp <- c()
c2_comp <- c()
for (idx in 1:n_rep) {
  print(idx)
  name_genes <- unique(df_wbc_tmp1$gene)  
  name_genes <- name_genes[sample(1:length(name_genes),
    length(name_genes), replace = TRUE)]
  
  df_wbc_tmp11 <- tibble(gene = name_genes) %>%
    left_join(df_wbc_tmp1,
      by = "gene")
  
  c1_comp <- c(c1_comp,
    (tidy(lm(s2 ~ term + pheno +study, df_wbc_tmp11)) %>%
        filter(grepl("term", term)))$estimate
  )
  
  name_genes <- unique(df_wbc_tmp2$gene)
  name_genes <- name_genes[sample(1:length(name_genes),
    length(name_genes), replace = TRUE)]
  
  df_wbc_tmp22 <- tibble(gene = name_genes) %>%
    left_join(df_wbc_tmp2,
      by = "gene")
  
  c2_comp <- c(c2_comp,
    (tidy(lm(s2 ~ term + pheno +study, df_wbc_tmp22)) %>%
        filter(grepl("term", term)))$estimate
  )
}

tmp_meta <- bind_rows(
  tidy(lm(s2 ~ term + pheno +study, df_wbc_tmp1)) %>%
    filter(grepl("term", term)),
  tidy(lm(s2 ~ term + pheno +study, df_wbc_tmp2)) %>%
    filter(grepl("term", term))
) 



tmp_meta$bs_se <- c(sd(c1_comp), sd(c2_comp))

tmp_meta %>%
  mutate(newp = pnorm(abs(estimate/bs_se), lower.tail = FALSE))

tmp_meta %>%
  mutate(weight = 1/(bs_se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cov_twas %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

mega_twas %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

mesusie_twas %>%
  filter(p.value < 0.05/23000) %>%
  summarize(n = n())

195-168
195-143
195/168
195/143
195 - (168+143)/2
195 /((168+143)/2)
195-168
195-143


# qqplot
qq_twas2 <-  df_twas %>%
  filter(term %in% c("sushie", "mega", "mesusie")) %>%
  filter(pheno %in% c("sex", "smoke")) %>%
  mutate(term = factor(term,
    levels = c("sushie", "mega", "mesusie"),
    labels = c("SuShiE", "SuSiE", "MESuSiE"))) %>%
  mutate(study = factor(study, 
    levels = c("mesa.mrna", "mesa.proteins", "genoa.mrna"),
    labels = c("TOPMed-MESA mRNA", "TOPMed-MESA Proteins", "GENOA mRNA")))

ggplot(qq_twas2, aes(sample = -log(p.value))) +
  stat_qq(aes(color = term), distribution = qexp) +
  stat_qq_line(distribution = qexp) +
  scale_color_manual(values = twas_colors) +
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

# ggsave("./plots/s28.png", width = p_width-1, height = p_height+3)

df_eds <- read_tsv(glue("{constraint_folder}/df_eds.tsv"))
df_pli <- read_tsv(glue("{constraint_folder}/df_pli_new.tsv"))
df_rvis <- read_tsv(glue("{constraint_folder}/df_rvis.tsv"))
df_shet <- read_tsv(glue("{constraint_folder}/df_shet_new.tsv"))

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

df_wk3 <- mesusie_twas %>%
  select(type = term, statistic, gene, pheno, study) %>%
  mutate(chisq_twas = statistic^2,
    trait = gene,
    gene = gsub("_.*", "", gene)) %>%
  left_join(df_scores %>%
      rename(gene = trait)) %>%
  filter(!is.na(value))

df_wk <- bind_rows(df_wk1, df_wk2, df_wk3)

df_asso <- df_wk %>%
  nest_by(type, score) %>%
  mutate(mod = list(lm(value ~ chisq_twas +pheno +study, data = data))) %>%
  reframe(tidy(mod)) %>%
  filter(term == "chisq_twas")  %>%
  mutate(score = factor(score, levels = c("pLI", "LOEUF", "s_het", "RVIS", "EDS")),
    type = factor(type, levels = c("sushie", "mega", "mesusie"))) %>%
  arrange(type, score)

pnorm(abs((sum(df_asso$statistic[1:5]^2) - sum(df_asso$statistic[6:10]^2))/sqrt(20)),
  lower.tail = FALSE)

pnorm(abs((sum(df_asso$statistic[1:5]^2) - sum(df_asso$statistic[11:15]^2))/sqrt(20)),
  lower.tail = FALSE)

# write_tsv(select(df_asso, estimate, statistic, p.value), "./tables/s15.tsv")

# validation
# lu et al
df_ref <- read_tsv(glue("{metadata_folder}/gencode.v34.gene.only.tsv.gz"))

load(glue("{validation_folder}/lu/twas.RData"))

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
lfs <- list.files(glue("{validation_folder}/kachuri/"),
  pattern="Result", full.names = TRUE)
kachuri_twas_sig <- lfs %>% map_df(read_csv) %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "WBC") %>%
  filter(abs(zscore) > sig_thresh) %>%
  distinct(gene, pheno)

mon_twas <- read_csv(glue("{validation_folder}/kachuri/MONO_AA_whole_blood.csv")) %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "MON") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

neu_twas <- read_csv(glue("{validation_folder}/kachuri/NEU_AA_whole_blood.csv")) %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "NEU") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

wbc_twas <- read_csv(glue("{validation_folder}/kachuri/WBC_AA_whole_blood.csv")) %>%
  mutate(gene = gsub("\\..+", "", gene),
    pheno = "WBC") %>%
  filter(abs(zscore) > sig_thresh) %>%
  select(gene, pheno)

kachuri_twas_sig2 <- bind_rows(mon_twas, neu_twas, wbc_twas)

# tapia et al 
library(readxl)
tapia_twas <- read_excel(glue("{validation_folder}/tapia/gepi22436-sup-0003-supptable2_twassignif.xlsx"), sheet = 1, skip = 3)

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

rowland_twas <- read_excel(glue("{validation_folder}/rowland/supplemental_tables_ddac011.xlsx"),sheet = 2) %>%
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
wen_twas_sig <- read_excel(glue("{validation_folder}/wen/Supplementary_Tables_finalSUB_R2.xlsx"), sheet = 4, skip = 1) %>%
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

mesusie_twas %>%
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

prop.test(c(110, 98), c(189, 162))
prop.test(c(110, 87), c(189, 141))

# main figure 5
sig_thresh <- qnorm(0.05/23000, lower.tail = FALSE)

main_twas1 <- cov_twas %>%
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
      "SuSiE", "MESuSiE",
      "Neither")))

all_levels <- c("Both",
  "SuShiE",
  "SuSiE", "MESuSiE",
  "Neither")

p1 <- ggplot(main_twas1, aes(x = susie, y = sushie, color = type)) +
  geom_point() +
  geom_abline(slope=1, intercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(values =
      c("#c994c7", "#1b9e77", "#e7298a", "#e6ab02", "lightgrey"),
    drop=FALSE) +
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

main_twas2 <- cov_twas %>%
  select(study, gene, pheno, sushie = statistic) %>%
  inner_join(mesusie_twas %>%
      select(study, gene, pheno, mesusie = statistic)) %>%
  mutate(type = ifelse(abs(sushie) > sig_thresh &
      abs(mesusie) > sig_thresh, "Both",
    ifelse(abs(sushie) > sig_thresh, "SuShiE",
      ifelse(abs(mesusie) > sig_thresh, "MESuSiE",
        "Neither"))),
    type = factor(type, levels = c("Both",
      "SuShiE",
      "MESuSiE",
      "Neither")))

p2 <- ggplot(main_twas2, aes(x = mesusie, y = sushie, color = type)) +
  geom_point() +
  geom_abline(slope=1, intercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#c994c7", "#1b9e77", "#e6ab02", "lightgrey")) +
  xlab(expression(bold("MESuSiE T/PWAS" ~ t))) +
  ylab(expression(bold("SuShiE T/PWAS" ~ t))) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    axis.title=element_text(face="bold", size = 7),
    axis.text = element_text(face = "bold", size = 6)) 


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

tt3 <- df_wk3 %>%
  group_by(Cate, score) %>%
  summarize(mvalue = mean(chisq_twas),
    se = sd(chisq_twas)/sqrt(n()),
    low_bound = mvalue - 1.96 *se,
    upp_bound = mvalue+1.96*se) %>%
  mutate(method = "MESuSiE")

tt <- bind_rows(tt1, tt2, tt3) %>%
  mutate(method = factor(method, levels = c("SuShiE", "SuSiE", "MESuSiE")),
    score = factor(score, levels = c("pLI", "LOEUF", "s_het", "RVIS", "EDS"),
      labels = c("bold(pLI)", "bold(LOEUF)", "bold(s[het])", "bold(RVIS)",
        "bold(EDS)")),
    Cate = factor(Cate, levels = c("Low", "Middle", "High")))

p3 <- ggplot(tt, aes(x = Cate, y = mvalue, fill = method, color = method)) +
  geom_point(position=position_dodge(width=0.5), shape=21, size=2) +
  geom_errorbar(aes(ymin = low_bound, ymax = upp_bound),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_fill_manual(values = c("#1b9e77", "#e7298a", "#e6ab02")) +
  scale_color_manual(values = c("#1b9e77", "#e7298a", "#e6ab02")) +
  ylab(expression(bold("Average T/PWAS "~ Chi^2))) +
  xlab("Constraint Group") +
  facet_grid(cols = vars(score), labeller = label_parsed) +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.position = "none",
    strip.text.x = element_text(face = "bold"),
    # legend.box.background = element_rect(colour = "black"),
    axis.title=element_text(face="bold", size = 7),
    axis.text = element_text(face = "bold", size = 6)) 

library(patchwork)
p1_legend <- cowplot::get_legend(p1)
try_design <- "
AB
AB
AB
AB
CC
CC
CC
CC
DD
"
p1+p2 + p3+p1_legend+ plot_layout(design = try_design, guides = "collect") &
  theme(legend.position="none")

# ggsave("./plots/p5.png", width = p_width-2, height = p_height+2)

tidy(lm(sushie ~susie + study, main_twas1))

tmp1 <- main_twas1 %>%
  mutate(sushie2 = sushie^2,
    susie2 = susie^2)

tidy(t.test(tmp1$sushie2, tmp1$susie2, alternative = "greater"))

tmp2 <- main_twas2 %>%
  mutate(sushie2 = sushie^2,
    mesusie2 = mesusie^2)

tidy(t.test(tmp2$sushie2, tmp2$mesusie2, alternative = "greater"))

