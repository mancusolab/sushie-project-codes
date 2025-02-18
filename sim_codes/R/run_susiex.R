library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

prefix <- args[1]
sim <- args[2]
locus <- args[3]
causal_file <- args[4]
out <- args[5]
L2 <- as.numeric(args[6])

df_ref <- read_tsv(causal_file)

path_cs <- glue("{prefix}/susiex.in.sim{sim}.locus{locus}.cs")
path_snp <- glue("{prefix}/susiex.in.sim{sim}.locus{locus}.snp")
if (file.exists(path_cs)) {
  df_cs2 <- read_tsv(path_cs, col_names = FALSE)
  if (df_cs2$X1[1] == "FAIL") {
    df_sens <- tibble("method" = "susiex",
      "sim" = sim,
      "locus" = locus,
      sens = 0)
    write_tsv(df_sens, glue("{out}.in.sim{sim}.locus{locus}.sens.tsv"))
    q()
  }
  
  df_sens <- tibble("method" = "susiex",
    "sim" = sim,
    "locus" = locus,
    sens = 1)
  write_tsv(df_sens, glue("{out}.in.sim{sim}.locus{locus}.sens.tsv"))
  
  if (df_cs2$X1[1] == "NULL") {
    df_res <- tibble(
      CSIndex = 1:L2,
      susiex = NA)
    
    df_cs_SNP <- c()
  } else {
    df_cs <- read_tsv(path_cs)
    df_res <- tibble(CSIndex = 1:L2) %>%
      left_join(df_cs %>%
          rename(CSIndex = CS_ID) %>%
          group_by(CSIndex) %>%
          summarize(susiex = n()),
        by = "CSIndex")
    
    df_cs_SNP <- df_cs$SNP
  }
  
  df_res$sim = sim
  df_res$locus = locus
  write_tsv(df_res, glue("{out}.in.sim{sim}.locus{locus}.cs.tsv"))
  
  if (file.exists(path_snp)){
    df_snp <- read_tsv(path_snp)
    n_cs <- sum(grepl("PIP", colnames(df_snp)))
    if (n_cs < L2) {
      for (idx in 1:(L2-n_cs)) {
        df_snp[[glue("PIP{idx}")]] <- rep(1/nrow(df_snp), nrow(df_snp))
      }
    }
    df_mat <- df_snp %>%
      select(contains("PIP")) %>%
      mutate(across(everything(), ~ 1 - .)) %>%
      rowwise() %>%
      mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
      ungroup()
    
    df_snp$pip_all <- df_mat$result 
    df_pip <- df_ref %>%
      left_join(df_snp %>% select(SNP, pip = pip_all), by = "SNP") %>%
      mutate(sim = sim,
        locus = locus,
        method = "susiex")
    df_pip$cali = as.numeric(df_ref$SNP %in% df_cs_SNP)

    write_tsv(df_pip, glue("{out}.in.sim{sim}.locus{locus}.pip.tsv"))
  }
  
}

