library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_2pop/sim3_locus233/other.in.sim3.locus233.rdata"

f_rdata <- args[1]
out <- args[2]
out2 <- args[3]
load(f_rdata)

# run MESuSiE
mesusie_res <- tryCatch({
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=L2,
    ancestry_weight = c(0, 0, 1),
    estimate_residual_variance =TRUE, max_iter =500)
  mesusie_res
}, error = function(e) {
  print("mesusie doesn't converge...")
  df_sens <- tibble("mesusie" = 0)
  df_sens$sim <- sim
  df_sens$locus <- locus
  write_tsv(df_sens, glue("{out}.sens.tsv"))
  q()
})


df_sens <- tibble("mesusie" = 1)
df_sens$sim <- sim
df_sens$locus <- locus
write_tsv(df_sens, glue("{out}.sens.tsv"))

df_tmp <- tibble(snp_index = 1:nrow(df_ld1))
for (idx in 1:L2) {
  df_tmp[[glue("PIP{idx}")]] <- rowSums(mesusie_res$alpha[[idx]])
}

df_mat <- df_tmp %>%
  select(contains("PIP")) %>%
  mutate(across(everything(), ~ 1 - .)) %>%
  rowwise() %>%
  mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
  ungroup()

df_causal$mesusie <- df_mat$result[df_causal$`SNPIndex_1based`]
df_causal$sim <- sim
df_causal$locus <- locus
df_causal$cali <- as.numeric(df_causal$`SNPIndex_1based` %in% unlist(mesusie_res$cs$cs))
write_tsv(df_causal, glue("{out}.pip.tsv"))


df_cs <- tibble("CSIndex" = seq(L2))
df_cs$sim <- sim
df_cs$locus <- locus
df_cs$mesusie <- NA

if (length(mesusie_res$cs$cs) != 0) {
  # for every credible set
  for (jdx in 1:length(mesusie_res$cs$cs)) {
    df_cs$mesusie[jdx] <- length(mesusie_res$cs$cs[[jdx]])
  }
}

write_tsv(df_cs, glue("{out}.cs.tsv"))


## additional
df_tmp <- tibble(snp_index = 1:nrow(df_ld1))
for (idx in 1:L2) {
  df_tmp[[glue("PIP{idx}")]] <- rowSums(mesusie_res$alpha[[idx]])
}

df_mat <- df_tmp %>%
  select(contains("PIP")) %>%
  mutate(across(everything(), ~ 1 - .)) %>%
  rowwise() %>%
  mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
  ungroup()

df_pip <- df_mat %>%
  mutate(SNPIndex_1based = row_number(),
    sim = sim,
    locus = locus) %>%
  select(sim, locus, SNPIndex_1based, mesusie = result)

write_tsv(df_pip, glue("{out2}.pip.tsv"))

if (length(mesusie_res$cs$cs) != 0) {
  # for every credible set
  df_cs <- tibble()
  for (jdx in 1:length(mesusie_res$cs$cs)) {
    df_cs <- df_cs %>%
      bind_rows(
        tibble(CSIndex = jdx, 
          SNPIndex_1based = mesusie_res$cs$cs[[jdx]]
        )
      )
  }
  df_cs$sim <- sim
  df_cs$locus <- locus
  
  write_tsv(df_cs, glue("{out2}.cs.tsv"))
}


# rho_res <- c()
# for (idx in 1:L2) {
#   res1 <- cbind(mesusie_res$mu1[[idx]][[1]], mesusie_res$mu1[[idx]][[2]])
#   res2 <- mesusie_res$mu1[[idx]][[3]]
#   res3 <- cbind(mesusie_res$mu1[[idx]][[1]] + mesusie_res$mu1[[idx]][[3]][,1],
#     mesusie_res$mu1[[idx]][[2]] + mesusie_res$mu1[[idx]][[3]][,2])
#   rho_res <- rho_res %>%
#     bind_rows(
#       tibble(value = as.vector(t(res1) %*% res1), num_orer = 1:4) %>%
#     mutate(CSIndex = idx,
#       type = "method1"),
#       tibble(value = as.vector(t(res2) %*% res2), num_orer = 1:4) %>%
#         mutate(CSIndex = idx,
#           type = "method2"),
#       tibble(value = as.vector(t(res3) %*% res3), num_orer = 1:4) %>%
#         mutate(CSIndex = idx,
#           type = "method3"),
#     )
# }

df_rho <- tibble("locus" = locus, "sim" = sim, CSIndex = 1:L2)
df_rho$est_rho <- NA
df_rho$frob <- NA
for (idx in 1:L2) {
  tmp_matrix <- mesusie_res$V[[idx]]
  tmp_rho <- tmp_matrix[1,2] / (sqrt(tmp_matrix[1,1]) * sqrt(tmp_matrix[2,2]))
  tmp_frob <- norm(tmp_matrix, type = "F")
  df_rho$est_rho[idx] <- tmp_rho
  df_rho$frob[idx] <- tmp_frob
}

df_rho <- df_rho %>%
  arrange(desc(frob)) %>%
  mutate(CSIndex = row_number())

write_tsv(df_rho, glue("{out2}.rho.tsv"))


