library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

f_rdata <- args[1]
out <- args[2]
load(f_rdata)


# run MESuSiE
mesusie_res <- tryCatch({
  mesusie_res <- MESuSiE::meSuSie_core(ld_list,ss_list, L=L2,
    ancestry_weight = c(0, 0, 1),
    estimate_residual_variance =TRUE, max_iter =500)
  mesusie_res
}, error = function(e) {
  print("mesusie cannot converge")
  res <- tibble("weights1" = rep(0, nrow(df_ss2)), "weights2" = rep(0, nrow(df_ss2)))
  write_tsv(res, out, col_names = FALSE)
  q()
})

res1 <- NULL
res2 <- NULL
for (idx in 1:L2) {
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

write_tsv(res, out, col_names = FALSE)

