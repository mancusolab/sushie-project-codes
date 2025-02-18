library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_2pop/sim4_locus468/other.in.sim4.locus468.rdata"

f_rdata <- args[1]
out <- args[2]
load(f_rdata)

# run MESuSiE
mesusie_res <- tryCatch({
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=10,
    estimate_residual_variance =TRUE, max_iter =500)
  mesusie_res
}, error = function(e) {
  print("mesusie doesn't converge...")
  q()
})

# get pip and weights
df_cs <- tibble()
if (length(mesusie_res$cs$cs) != 0) {
  # for every credible set
  for (jdx in 1:length(mesusie_res$cs$cs)) {
    df_tmp <- tibble(SNPIndex = mesusie_res$cs$cs[[jdx]], CSIndex = jdx)
    df_cs <- bind_rows(df_cs, df_tmp)
  }
}


if (nrow(df_cs) != 0) {
  res1 <- NULL
  res2 <- NULL
  df_alpha <- NULL
  for (idx in 1:10) {
    
    if (idx %in% unique(df_cs$CSIndex)) {
      tmp_alpha <- tibble(x = rowSums(mesusie_res$alpha[[idx]])) 
      colnames(tmp_alpha) <- glue("alpha_l{idx}")
      df_alpha <- bind_cols(df_alpha, tmp_alpha)
    }
    
    for (jdx in 1:3) {
      if (jdx == 1) {
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] * mesusie_res$mu1[[idx]][[jdx]])
      } else if (jdx == 2){
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] * mesusie_res$mu1[[idx]][[jdx]])
      } else {
        res1 <- cbind(res1, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,1])
        res2 <- cbind(res2, mesusie_res$alpha[[idx]][,jdx] *
            mesusie_res$mu1[[idx]][[jdx]][,2])
      }
    }
  }
  
  df_tmp <- tibble(snp_index = 1:nrow(df_ref))
  for (idx in 1:10) {
    df_tmp[[glue("PIP{idx}")]] <- rowSums(mesusie_res$alpha[[idx]])
  }
  df_mat <- df_tmp %>%
    select(contains("PIP")) %>%
    mutate(across(everything(), ~ 1 - .)) %>%
    rowwise() %>%
    mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
    ungroup()
  
  df_res <- df_cs %>%
    left_join(df_ref %>%
        select(chrom, snp, pos, a0, a1) %>%
        mutate(trait = trait_name)%>%
        rownames_to_column() %>%
        rename(SNPIndex = rowname) %>%
        mutate(SNPIndex = as.numeric(SNPIndex)),
      by = "SNPIndex") %>%
    mutate(method = "mesusie")
  
  alpha_list <- c()
  for (idx in unique(df_res$CSIndex)) {
    n_tmp <- df_res %>%
      filter(CSIndex == idx) %>%
      select(SNPIndex) %>%
      unlist() %>%
      as.numeric()
    alpha_list <- c(alpha_list, df_alpha[[glue("alpha_l{idx}")]][n_tmp])
  }
  df_res$alpha <- alpha_list
  df_res$pip_all <- df_mat$result[df_res$SNPIndex]
  
  write_tsv(df_res, glue("{out}/cs/mesusie.{trait_name}.cs.tsv"))
  
  res <- tibble(as.data.frame(cbind(rowSums(res1), rowSums(res2))))
  df_ref$pop1_weights <- rowSums(res1)
  df_ref$pop2_weights <- rowSums(res2)
  df_ref$trait <- trait_name
  df_ref$method <- "mesusie"
  df_ref <- bind_cols(df_ref, df_alpha)
  
  df_ref$pip_all <- df_mat$result
  
  write_tsv(df_ref, glue("{out}/weights/mesusie.{trait_name}.weights.tsv"))
}