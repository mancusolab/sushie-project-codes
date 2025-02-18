library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_2pop/sim4_locus468/other.in.sim4.locus468.rdata"

f_rdata <- args[1]
out <- args[2]
load(f_rdata)

# run MESuSiE
mesusie_res <- tryCatch({
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=L2, estimate_residual_variance =TRUE, max_iter =500)
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



