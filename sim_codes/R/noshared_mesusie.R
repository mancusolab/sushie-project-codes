library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

# f_rdata <- "/scratch1/zeyunlu/tmp_noshared/sim4_locus497/other.in.sim4.locus497.rdata"

f_rdata <- args[1]
out <- args[2]
load(f_rdata)

# run MESuSiE
mesusie_res <- tryCatch({
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=L2,
    ancestry_weight = c(1, 1, 0),
    estimate_residual_variance =TRUE, max_iter =500)
  mesusie_res
}, error = function(e) {
  print("mesusie doesn't converge...")
  q()
})


df_fdr <- tibble("method" = "mesusie")

df_cs <- tibble()
if (length(mesusie_res$cs$cs) != 0) {
  # for every credible set
  for (jdx in 1:length(mesusie_res$cs$cs)) {
      df_tmp <- tibble(SNPIndex = mesusie_res$cs$cs[[jdx]], CSIndex = jdx)
      df_cs <- bind_rows(df_cs, df_tmp)
  }
}

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

if (nrow(df_cs)!=0) {
  tmp_cs <- df_cs %>%
    filter(SNPIndex %in% df_causal$`SNPIndex_1based`) %>%
    distinct(CSIndex)
  cs_num <- length(unique(df_cs$CSIndex))
  as_in_num <- nrow(tmp_cs)
  fdr_high_num <- sum(df_mat$result > 0.95)
  as_high_num <- sum(df_mat$result > 0.95 &
      1:length(df_mat$result) %in% df_causal$`SNPIndex_1based`)
} else {
  cs_num <- 0
  as_in_num <- 0
  fdr_high_num <- 0
  as_high_num <- 0
}
df_fdr$fdr_cs <- cs_num
df_fdr$as_in <- as_in_num
df_fdr$fdr_high <- fdr_high_num
df_fdr$as_high <- as_high_num
df_fdr$sim <- sim
df_fdr$locus <- locus

write_tsv(df_fdr, glue("{out}.fdr.tsv"))


df_causal$pip <- df_mat$result[df_causal$`SNPIndex_1based`]
df_causal$sim <- sim
df_causal$locus <- locus
df_causal$method <- "mesusie"
df_causal$ancestry <- c(rep("1", L2/2), rep("2", L2/2))
df_causal$cali <- as.numeric(df_causal$`SNPIndex_1based` %in% mesusie_res$cs$cs)
write_tsv(df_causal, glue("{out}.pip.tsv"))
