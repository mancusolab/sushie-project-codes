library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

path_cs <- args[1]
path_snp <- args[2]
trait_name <- args[3]
chrom_name <- args[4]
out <- args[5]

if (file.exists(path_cs)) {
  
  df_cs2 <- read_tsv(path_cs, col_names = FALSE)
  
  if (df_cs2$X1[1] != "FAIL" & df_cs2$X1[1] != "NULL") {
    df_cs <- read_tsv(path_cs)
    df_snp <- read_tsv(path_snp)
    
    n_cs <- length(unique(df_cs$CS_ID))
    df_snp_bp <- df_snp
    if (n_cs < 10) {
      for (idx in 1:(10-n_cs)) {
        df_snp_bp[[glue("PIP{idx}")]] <- rep(1/nrow(df_snp_bp), nrow(df_snp_bp))
      }
    }
    
    df_mat <- df_snp_bp %>%
      select(contains("PIP")) %>%
      mutate(across(everything(), ~ 1 - .)) %>%
      rowwise() %>%
      mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
      ungroup()
    
    df_mat2 <- df_snp %>%
      select(SNP)
    df_mat2$pip_all <- df_mat$result
    
    df_cs2 <- df_cs %>%
      mutate(trait = trait_name) %>%
      left_join(df_mat2, by = "SNP") %>%
      select(CSIndex = CS_ID, snp = SNP, pos = BP, alpha = CS_PIP,
        pip_all, trait) %>%
      mutate(method = "susiex")
    
    write_tsv(df_cs2, glue("{out}/cs/susiex.{trait_name}.cs.tsv"))
    
    df_res <- df_snp %>%
      select(SNP, contains("PIP")) %>%
      left_join(df_mat2, by = "SNP") %>%
      mutate(trait = trait_name,
        chrom = chrom_name,
        method = "susiex")
    
    pip_columns <- grep("^PIP\\(CS\\d+\\)$", names(df_res))
    new_names <- paste0("alpha_l", seq_along(pip_columns))
    names(df_res)[pip_columns] <- new_names
    
    write_tsv(df_res, glue("{out}/weights/susiex.{trait_name}.weights.tsv"))
  }
}
