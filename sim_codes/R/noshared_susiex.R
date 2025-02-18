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
    q()
  }
  
  if (df_cs2$X1[1] == "NULL") {
    df_res <- tibble("method" = "susiex",
      "sim" = sim,
      "locus" = locus,
      fdr_cs = 0,
      as_in = 0,
      fdr_high = 0,
      as_high = 0)
    
  } else {
    df_cs <- read_tsv(path_cs)
    
    if (file.exists(path_snp)) {
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
      df_res <- df_ref %>%
        left_join(df_snp %>% select(SNP, pip = pip_all), by = "SNP") %>%
        mutate(sim = sim,
          locus = locus,
          method = "susiex")
      df_res$cali = as.numeric(df_ref$SNP %in% df_cs$SNP)
      df_res$ancestry <- c(rep("1", L2/2), rep("2",  L2/2))
      
      write_tsv(df_res, glue("{out}.in.sim{sim}.locus{locus}.pip.tsv"))
    }
    
    tmp_cs <- df_cs %>%
      filter(SNP %in% df_ref$SNP) %>%
      distinct(CS_ID)
    
    df_cs <- df_cs %>%
      left_join(df_snp %>% select(SNP, pip_all), by = "SNP")
    
    df_res <- tibble("method" = "susiex",
      "sim" = sim,
      "locus" = locus,
      fdr_cs = length(unique(df_cs$CS_ID)),
      as_in = nrow(tmp_cs),
      fdr_high = sum(df_cs$pip_all > 0.95),
      as_high = sum(df_cs$pip_all > 0.95 & df_cs$SNP %in% df_ref$SNP))
    
  }
  
  write_tsv(df_res, glue("{out}.in.sim{sim}.locus{locus}.fdr.tsv"))
  
}

