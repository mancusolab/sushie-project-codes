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
res1 <- NULL
res2 <- NULL
for (idx in 1:10) {
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

res <- tibble(as.data.frame(cbind(rowSums(res1), rowSums(res2))))
df_ref$pop1_weights <- rowSums(res1)
df_ref$pop2_weights <- rowSums(res2)
df_ref$trait <- trait_name
df_ref$method <- "mesusie"
write_tsv(df_ref, out)

