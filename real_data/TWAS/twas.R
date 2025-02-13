library(tidyverse)
library(broom)
library(glue)

perform_rint <- function(pheno) {
  dt_rank <- rank(pheno)
  dt_n <- length(pheno)
  dt_irn <- qnorm(dt_rank / (dt_n + 1))
  return(dt_irn)
}

args <- commandArgs(trailingOnly = TRUE)

chrom_name <- as.numeric(args[1])
gene_name <- args[2]
study_name <- args[3]
ge_file <- args[4]
pheno_prefix <- args[5]
pt_file <- args[6]
exc_file <- args[7]
out_file <- args[8]

print(gene_name)

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

traits <- c("BAS", "EOS", "LYM", "MON", "NEU", "WBC")

final_res <- tibble()
for (trait in traits) {
  print(trait)
  tmp_df_pheno <- read_table(glue("{pheno_prefix}/{trait}.tsv"),
    col_names=FALSE, col_types = cols()) %>%
    # make sex binary 0 and 1
    mutate(X4 = X4 -1)
  
  colnames(tmp_df_pheno) <-
    c("FID", "IID", "PHENO", "SEX", "AGE", "AGE2", paste0("PC", 1:10))
  
  df_wk <- tmp_df_pheno %>%
    inner_join(eur_pt, by = "IID") %>%
    left_join(df_ge %>%
        select(IID,
          sushie = sushie_pop1_AVG,
          indep = indep_pop1_AVG,
          meta = meta_pop1_AVG,
          mega = mega_weight_AVG,
          mesusie = mesusie_pop1_AVG,
          enet = enet_weight_AVG,
          lasso = lasso_weight_AVG,
          gblup = gblup_weight_AVG) %>%
        mutate(sushie = scale(sushie),
          indep = scale(indep),
          meta = scale(meta),
          mega = scale(mega),
          mesusie = scale(mesusie),
          enet = scale(enet),
          lasso = scale(lasso),
          gblup = scale(gblup)),
      by = "IID") %>%
    bind_rows(tmp_df_pheno %>%
        inner_join(afr_pt, by = "IID") %>%
        left_join(df_ge %>%
            select(IID,
              sushie = sushie_pop2_AVG,
              indep = indep_pop2_AVG,
              meta = meta_pop2_AVG,
              mega = mega_weight_AVG,
              mesusie = mesusie_pop2_AVG,
              enet = enet_weight_AVG,
              lasso = lasso_weight_AVG,
              gblup = gblup_weight_AVG) %>%
            mutate(sushie = scale(sushie),
              indep = scale(indep),
              meta = scale(meta),
              mega = scale(mega),
              mesusie = scale(mesusie),
              enet = scale(enet),
              lasso = scale(lasso),
              gblup = scale(gblup)),
          by = "IID"))
  
  if (study_name != "genoa.mrna") {
    df_wk <- df_wk %>%
      bind_rows(tmp_df_pheno %>%
          inner_join(his_pt, by = "IID") %>%
          left_join(df_ge %>%
              select(IID,
                sushie = sushie_pop3_AVG,
                indep = indep_pop3_AVG,
                meta = meta_pop3_AVG,
                mega = mega_weight_AVG,
                mesusie = mesusie_pop3_AVG,
                enet = enet_weight_AVG,
                lasso = lasso_weight_AVG,
                gblup = gblup_weight_AVG) %>%
              mutate(sushie = scale(sushie),
                indep = scale(indep),
                meta = scale(meta),
                mega = scale(mega),
                mesusie = scale(mesusie),
                enet = scale(enet),
                lasso = scale(lasso),
                gblup = scale(gblup)),
            by = "IID"))
  }
  
  pheno_mean <- mean(df_wk$PHENO, na.rm = TRUE)
  pheno_std <- sd(df_wk$PHENO, na.rm = TRUE)
  
  df_wk <- df_wk %>%
    filter(PHENO >= pheno_mean - 3 * pheno_std & PHENO <= pheno_mean + 3 * pheno_std)
  
  df_res <- lm(PHENO ~ SEX + AGE + AGE2 + PC1 + PC2 + PC3 + PC4
    + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk)$residuals
  rint_res <- perform_rint(df_res)
  
  tmp_all_res <- tibble()
  for (idx in c("sushie", "indep", "meta", "mega", "mesusie", "enet", "lasso", "gblup")) {
    df_tmp <- df_wk[c(idx, "SEX", "AGE", "AGE2", paste0("PC", 1:10))]
    colnames(df_tmp) <- c("method", "SEX", "AGE", "AGE2", paste0("PC", 1:10))
    if (sum(df_tmp$method, na.rm = TRUE) == 0 |
        all(is.na(df_tmp$method))) { 
      next
    }
    
    if (any(is.na(df_tmp$method))) {
      tmp_rint_res <- rint_res[!is.na(df_tmp$method)]
      df_tmp <- df_tmp %>%
        filter(!is.na(method))
    } else {
      tmp_rint_res <- rint_res
    }
    
    method_res <- lm(method ~ ., data = df_tmp)$residuals
    res_mod <- tidy(lm(tmp_rint_res ~ method_res)) %>%
      filter(term != "(Intercept)") %>%
      mutate(term = idx)
    tmp_all_res <- tmp_all_res %>%
      bind_rows(res_mod)
  }
  
  num_eur <- sum(df_wk$IID %in% eur_pt$IID)
  num_afr <- sum(df_wk$IID %in% afr_pt$IID)
  
  if (study_name != "genoa.mrna") {
    num_his <- sum(df_wk$IID %in% his_pt$IID)  
  } else {
    num_his <- 0
  }
  
  final_res <- final_res %>%
    bind_rows(
      tmp_all_res %>%
        mutate(gene = gene_name,
          pheno = trait,
          study = study_name,
          n_eur = num_eur,
          n_afr = num_afr,
          n_his = num_his))
}

# smoke TWAS
print("smoke TWAS")
tmp_df_pheno <- read_table(glue("{pheno_prefix}/smoke.tsv"),
  col_names=FALSE, col_types = cols()) %>%
  mutate(X4 = X4 -1)

colnames(tmp_df_pheno) <-
  c("IID","FID", "PHENO", "SEX", "AGE", "AGE2", paste0("PC", 1:10))

df_wk <- tmp_df_pheno %>%
  inner_join(eur_pt, by = "IID") %>%
  left_join(df_ge %>%
      select(IID,
        sushie = sushie_pop1_AVG,
        indep = indep_pop1_AVG,
        meta = meta_pop1_AVG,
        mega = mega_weight_AVG,
        mesusie = mesusie_pop1_AVG,
        enet = enet_weight_AVG,
        lasso = lasso_weight_AVG,
        gblup = gblup_weight_AVG) %>%
      mutate(sushie = scale(sushie),
        indep = scale(indep),
        meta = scale(meta),
        mega = scale(mega),
        mesusie = scale(mesusie),
        enet = scale(enet),
        lasso = scale(lasso),
        gblup = scale(gblup)),
    by = "IID") %>%
  bind_rows(tmp_df_pheno %>%
      inner_join(afr_pt, by = "IID") %>%
      left_join(df_ge %>%
          select(IID,
            sushie = sushie_pop2_AVG,
            indep = indep_pop2_AVG,
            meta = meta_pop2_AVG,
            mega = mega_weight_AVG,
            mesusie = mesusie_pop2_AVG,
            enet = enet_weight_AVG,
            lasso = lasso_weight_AVG,
            gblup = gblup_weight_AVG) %>%
          mutate(sushie = scale(sushie),
            indep = scale(indep),
            meta = scale(meta),
            mega = scale(mega),
            mesusie = scale(mesusie),
            enet = scale(enet),
            lasso = scale(lasso),
            gblup = scale(gblup)),
        by = "IID"))

if (study_name != "genoa.mrna") {
  df_wk <- df_wk %>%
    bind_rows(tmp_df_pheno %>%
        inner_join(his_pt, by = "IID") %>%
        left_join(df_ge %>%
            select(IID,
              sushie = sushie_pop3_AVG,
              indep = indep_pop3_AVG,
              meta = meta_pop3_AVG,
              mega = mega_weight_AVG,
              mesusie = mesusie_pop3_AVG,
              enet = enet_weight_AVG,
              lasso = lasso_weight_AVG,
              gblup = gblup_weight_AVG) %>%
            mutate(sushie = scale(sushie),
              indep = scale(indep),
              meta = scale(meta),
              mega = scale(mega),
              mesusie = scale(mesusie),
              enet = scale(enet),
              lasso = scale(lasso),
              gblup = scale(gblup)),
          by = "IID"))
}


df_res <- glm(PHENO ~ SEX + AGE + AGE2 + PC1 + PC2 + PC3 + PC4
  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk, family = binomial(link="logit"))$residuals
rint_res <- perform_rint(df_res)

tmp_all_res <- tibble()
for (idx in c("sushie", "indep", "meta", "mega", "mesusie", "enet", "lasso", "gblup")) {
  df_tmp <- df_wk[c(idx, "SEX", "AGE", "AGE2", paste0("PC", 1:10))]
  colnames(df_tmp) <- c("method", "SEX", "AGE", "AGE2", paste0("PC", 1:10))
  if (sum(df_tmp$method, na.rm = TRUE) == 0 |
      all(is.na(df_tmp$method))) { 
    next
  }
  
  if (any(is.na(df_tmp$method))) {
    tmp_rint_res <- rint_res[!is.na(df_tmp$method)]
    df_tmp <- df_tmp %>%
      filter(!is.na(method))
  } else {
    tmp_rint_res <- rint_res
  }
  
  method_res <- lm(method ~ ., data = df_tmp)$residuals
  res_mod <- tidy(lm(tmp_rint_res ~ method_res)) %>%
    filter(term != "(Intercept)") %>%
    mutate(term = idx)
  tmp_all_res <- tmp_all_res %>%
    bind_rows(res_mod)
}

num_eur <- sum(df_wk$IID %in% eur_pt$IID)
num_afr <- sum(df_wk$IID %in% afr_pt$IID)

if (study_name != "genoa.mrna") {
  num_his <- sum(df_wk$IID %in% his_pt$IID)  
} else {
  num_his <- 0
}

final_res <- final_res %>%
  bind_rows(
    tmp_all_res %>%
      mutate(gene = gene_name,
        pheno = "smoke",
        study = study_name,
        n_eur = num_eur,
        n_afr = num_afr,
        n_his = num_his))

# sex TWAS
print("SEX TWAS")
tmp_df_pheno <- read_table(glue("{pheno_prefix}/WBC.tsv"),
  col_names=FALSE, col_types = cols()) %>%
  mutate(X4 = X4 -1)

colnames(tmp_df_pheno) <- c("FID", "IID", "WBC", "PHENO", "AGE", "AGE2", paste0("PC", 1:10))

df_wk <- tmp_df_pheno %>%
  select(IID, PHENO, AGE, AGE2, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
  inner_join(eur_pt, by = "IID") %>%
  left_join(df_ge %>%
      select(IID,
        sushie = sushie_pop1_AVG,
        indep = indep_pop1_AVG,
        meta = meta_pop1_AVG,
        mega = mega_weight_AVG,
        mesusie = mesusie_pop1_AVG,
        enet = enet_weight_AVG,
        lasso = lasso_weight_AVG,
        gblup = gblup_weight_AVG) %>%
      mutate(sushie = scale(sushie),
        indep = scale(indep),
        meta = scale(meta),
        mega = scale(mega),
        mesusie = scale(mesusie),
        enet = scale(enet),
        lasso = scale(lasso),
        gblup = scale(gblup)),
    by = "IID") %>%
  bind_rows(tmp_df_pheno %>%
      select(IID, PHENO, AGE, AGE2, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
      inner_join(afr_pt, by = "IID") %>%
      left_join(df_ge %>%
          select(IID,
            sushie = sushie_pop2_AVG,
            indep = indep_pop2_AVG,
            meta = meta_pop2_AVG,
            mega = mega_weight_AVG,
            mesusie = mesusie_pop2_AVG,
            enet = enet_weight_AVG,
            lasso = lasso_weight_AVG,
            gblup = gblup_weight_AVG) %>%
          mutate(sushie = scale(sushie),
            indep = scale(indep),
            meta = scale(meta),
            mega = scale(mega),
            mesusie = scale(mesusie),
            enet = scale(enet),
            lasso = scale(lasso),
            gblup = scale(gblup)),
        by = "IID"))

if (study_name != "genoa.mrna") {
  df_wk <- df_wk %>%
    bind_rows(tmp_df_pheno %>%
        select(IID, PHENO, AGE, AGE2, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>%
        inner_join(his_pt, by = "IID") %>%
        left_join(df_ge %>%
            select(IID,
              sushie = sushie_pop3_AVG,
              indep = indep_pop3_AVG,
              meta = meta_pop3_AVG,
              mega = mega_weight_AVG,
              mesusie = mesusie_pop3_AVG,
              enet = enet_weight_AVG,
              lasso = lasso_weight_AVG,
              gblup = gblup_weight_AVG) %>%
            mutate(sushie = scale(sushie),
              indep = scale(indep),
              meta = scale(meta),
              mega = scale(mega),
              mesusie = scale(mesusie),
              enet = scale(enet),
              lasso = scale(lasso),
              gblup = scale(gblup)),
          by = "IID"))
}

df_res <- glm(PHENO ~ AGE + AGE2 + PC1 + PC2 + PC3 + PC4
  + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = df_wk,
  family = binomial(link="logit"))$residuals
rint_res <- perform_rint(df_res)

tmp_all_res <- tibble()
for (idx in c("sushie", "indep", "meta", "mega", "mesusie", "enet", "lasso", "gblup")) {
  df_tmp <- df_wk[c(idx, "AGE", "AGE2", paste0("PC", 1:10))]
  colnames(df_tmp) <- c("method", "AGE", "AGE2", paste0("PC", 1:10))
  if (sum(df_tmp$method, na.rm = TRUE) == 0 |
      all(is.na(df_tmp$method))) { 
    next
  }
  
  if (any(is.na(df_tmp$method))) {
    tmp_rint_res <- rint_res[!is.na(df_tmp$method)]
    df_tmp <- df_tmp %>%
      filter(!is.na(method))
  } else {
    tmp_rint_res <- rint_res
  }
  
  method_res <- lm(method ~ ., data = df_tmp)$residuals
  res_mod <- tidy(lm(tmp_rint_res ~ method_res)) %>%
    filter(term != "(Intercept)") %>%
    mutate(term = idx)
  tmp_all_res <- tmp_all_res %>%
    bind_rows(res_mod)
}

num_eur <- sum(df_wk$IID %in% eur_pt$IID)
num_afr <- sum(df_wk$IID %in% afr_pt$IID)

if (study_name != "genoa.mrna") {
  num_his <- sum(df_wk$IID %in% his_pt$IID)  
} else {
  num_his <- 0
}

final_res <- final_res %>%
  bind_rows(
    tmp_all_res %>%
      mutate(gene = gene_name,
        pheno = "sex",
        study = study_name,
        n_eur = num_eur,
        n_afr = num_afr,
        n_his = num_his))

write_tsv(final_res, out_file)
