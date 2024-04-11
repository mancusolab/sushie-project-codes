library(tidyverse)
library(broom)

perform_rint <- function(pheno) {
  dt_rank <- rank(pheno)
  dt_n <- length(pheno)
  dt_irn <- qnorm(dt_rank / (dt_n + 1))
  return(dt_irn)
}

args <- commandArgs(trailingOnly = TRUE)

chrom_name <- args[1]
gene_name <- args[2]
study_name <- args[3]
ge_file <- args[4]
pheno_prefix <- args[5]
pheno_sufix <- args[6]
pt_file <- args[7]
exc_file <- args[8]
out_file <- args[9]

# read in exc file
exc_pt <- read_tsv(exc_file, col_types = cols())
colnames(exc_pt) <- c("X1", "IID")

# read in ancestry information file
df_ans <- read_tsv(pt_file, col_types = cols())

eur_pt <- df_ans %>%
  filter(!research_id %in% exc_pt$IID) %>%
  filter(ancestry_pred == "eur") %>%
  select(IID = research_id)

afr_pt <- df_ans %>%
  filter(!research_id %in% exc_pt$IID) %>%
  filter(ancestry_pred == "afr")  %>%
  select(IID = research_id)

if (study_name != "genoa.mrna") {
  his_pt <- df_ans %>%
    filter(!research_id %in% exc_pt$IID) %>%
    filter(ancestry_pred == "amr")  %>%
    select(IID = research_id)
}

# read in gene expression file
df_ge <- read_tsv(ge_file, col_types = cols())

traits <- c("BAS", "EOS", "HCT", "HGB", "LYM",
  "MCH", "MCHC", "MCV", "MON", "MPV", "NEU",
  "PLT", "RBC", "RDW", "WBC")

final_res <- tibble()
for (trait in traits) {
  tmp_df_pheno <- read_tsv(paste0(pheno_prefix, trait, pheno_sufix), col_names=FALSE,
    col_types = cols())
  colnames(tmp_df_pheno) <- c("row", "FID", "IID", "PHENO", "SEX", "AGE", "AGE2", paste0("PC", 1:10))
  
  df_wk <- tmp_df_pheno %>%
    inner_join(eur_pt, by = "IID") %>%
    left_join(df_ge %>%
        select(IID, ancestry = EUR_weight_AVG, mega = mega_weight_AVG) %>%
        mutate(ancestry = scale(ancestry),
          mega = scale(mega)),
      by = "IID") %>%
    bind_rows(tmp_df_pheno %>%
        inner_join(afr_pt, by = "IID") %>%
        left_join(df_ge %>%
            select(IID, ancestry = AFR_weight_AVG, mega = mega_weight_AVG) %>%
            mutate(ancestry = scale(ancestry),
              mega = scale(mega)),
          by = "IID"))
  
  if (study_name != "genoa.mrna") {
    df_wk <- df_wk %>%
      bind_rows(tmp_df_pheno %>%
          inner_join(his_pt, by = "IID") %>%
          left_join(df_ge %>%
              select(IID, ancestry = HIS_weight_AVG, mega = mega_weight_AVG) %>%
              mutate(ancestry = scale(ancestry),
                mega = scale(mega)),
            by = "IID"))
  }
  
  pheno_mean <- mean(df_wk$PHENO, na.rm = TRUE)
  pheno_std <- sd(df_wk$PHENO, na.rm = TRUE)
  
  df_wk <- df_wk %>%
    filter(PHENO >= pheno_mean - 3 * pheno_std & PHENO <= pheno_mean + 3 * pheno_std)
  
  df_res <- lm(PHENO ~ SEX + AGE + AGE2 + PC1 + PC2 + PC3 + PC4
    + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk)$residuals
  
  df_anc <- lm(ancestry ~ SEX + AGE + AGE2 + PC1 + PC2 + PC3 + PC4
    + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk)$residuals
  
  df_meg <- lm(mega ~ SEX + AGE + AGE2 + PC1 + PC2 + PC3 + PC4
    + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk)$residuals
  
  rint_res <- perform_rint(df_res)
  
  res_mod1 <- tidy(lm(rint_res ~ df_anc)) %>%
    filter(term != "(Intercept)")
  
  res_mod2 <- tidy(lm(rint_res ~ df_meg)) %>%
    filter(term != "(Intercept)")
  
  num_eur <- sum(df_wk$IID %in% eur_pt$IID)
  num_afr <- sum(df_wk$IID %in% afr_pt$IID)
  
  if (study_name != "genoa.mrna") {
    num_his <- sum(df_wk$IID %in% his_pt$IID)  
  } else {
    num_his <- 0
  }
  
  final_res <- final_res %>%
    bind_rows(
      bind_rows(res_mod1, res_mod2) %>%
        mutate(gene = gene_name,
          pheno = trait,
          study = study_name,
          n_eur = num_eur,
          n_afr = num_afr,
          n_his = num_his))
}

write_tsv(final_res, out_file)


