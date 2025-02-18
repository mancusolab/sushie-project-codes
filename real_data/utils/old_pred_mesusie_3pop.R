library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_2pop/sim4_locus468/other.in.sim4.locus468.rdata"

tmp_path <- args[1]
trait_name <- args[2]
output_file <- args[3]

all_predy <- c()
all_realy <- c()
for (cv in 1:5) {
  tmp_ss1 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans0.ss.tsv"))
  tmp_ss2 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans1.ss.tsv"))
  tmp_ss3 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans2.ss.tsv"))
  df_ld1 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans0.ld.tsv"), col_names = FALSE) %>%
    as.matrix()
  df_ld2 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans1.ld.tsv"), col_names = FALSE) %>%
    as.matrix()
  df_ld3 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans2.ld.tsv"), col_names = FALSE) %>%
    as.matrix()
  
  df_ss1 <- tmp_ss1 %>%
    mutate(snp = SNP,
      SNP = paste0("SNP", POS)) %>%
    select(snp, SNP, Beta, Se, Z, N, POS) %>%
    column_to_rownames(var = "snp") %>%
    as.data.frame()
  
  df_ss2 <- tmp_ss2 %>%
    mutate(snp = SNP,
      SNP = paste0("SNP", POS)) %>%
    select(snp, SNP, Beta, Se, Z, N, POS) %>%
    column_to_rownames(var = "snp") %>%
    as.data.frame()
  
  df_ss3 <- tmp_ss3 %>%
    mutate(snp = SNP,
      SNP = paste0("SNP", POS)) %>%
    select(snp, SNP, Beta, Se, Z, N, POS) %>%
    column_to_rownames(var = "snp") %>%
    as.data.frame()
  
  ss_list <- list("pop1" = df_ss1, "pop2" = df_ss2, "pop3" = df_ss3)
  
  # it might lose precision when pandas output this
  new_matrix <- df_ld1
  diag(new_matrix) <- 1
  new_matrix[upper.tri(new_matrix)] <- df_ld1[upper.tri(df_ld1)]
  new_matrix[lower.tri(new_matrix)] <- t(df_ld1)[lower.tri(new_matrix)]
  df_ld1 <- new_matrix
  colnames(df_ld1) <- df_ss1$SNP
  rownames(df_ld1) <- rownames(df_ss1)
  
  new_matrix <- df_ld2
  diag(new_matrix) <- 1
  new_matrix[upper.tri(new_matrix)] <- df_ld2[upper.tri(df_ld2)]
  new_matrix[lower.tri(new_matrix)] <- t(df_ld2)[lower.tri(new_matrix)]
  df_ld2 <- new_matrix
  colnames(df_ld2) <- df_ss2$SNP
  rownames(df_ld2) <- rownames(df_ss2)
  
  new_matrix <- df_ld3
  diag(new_matrix) <- 1
  new_matrix[upper.tri(new_matrix)] <- df_ld3[upper.tri(df_ld3)]
  new_matrix[lower.tri(new_matrix)] <- t(df_ld3)[lower.tri(new_matrix)]
  df_ld3 <- new_matrix
  colnames(df_ld3) <- df_ss3$SNP
  rownames(df_ld3) <- rownames(df_ss3)
  
  ld_list <- list("pop1" = df_ld1, "pop2" = df_ld2, "pop3" = df_ld3)
  
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=10,
    estimate_residual_variance =TRUE, max_iter =500)
  
  res1 <- NULL
  res2 <- NULL
  res3 <- NULL
  for (idx in 1:10) {
    for (jdx in 1:7) {
      if (jdx == 1) {
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] * mesusie_res$mu1[[idx]][[jdx]])
      } else if (jdx == 2){
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] * mesusie_res$mu1[[idx]][[jdx]])
      } else if (jdx == 3){
        res3 <- cbind(res3, mesusie_res$alpha[[idx]][,jdx] * mesusie_res$mu1[[idx]][[jdx]])
      } else if (jdx == 4){
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,1])
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,2])
      } else if (jdx == 5){
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,1])
        res3 <- cbind(res3, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]])
      } else if (jdx == 6) {
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,2])
        res3 <- cbind(res3, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]])
      } else {
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,1])
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,2])
        res3 <- cbind(res3, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]])
      }
    }
  }
  for (kdx in 1:3) {
    if (kdx == 1) {
      res_weight <- res1
    } else if (kdx == 2) {
      res_weight <- res2
    } else {
      res_weight <- res3
    }
    
    tmp_x1 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans{kdx-1}.validx.tsv"),
      col_names = FALSE) %>%
      as.matrix()
    tmp_y1 <- read_tsv(glue("{tmp_path}/group{cv-1}.ans{kdx-1}.validy.tsv"),
      col_names = FALSE)
    
    all_predy <- c(all_predy, tmp_x1 %*% rowSums(res_weight))
    all_realy <- c(all_realy, tmp_y1$X1)
  }
  
}

final_res <- tibble("mesusie" = summary(lm(all_predy ~ all_realy))$r.squared,
  trait = trait_name)

write_tsv(final_res, output_file)

