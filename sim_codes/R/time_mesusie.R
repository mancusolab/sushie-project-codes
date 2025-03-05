library(tidyverse)
library(glue)
library(microbenchmark)
args <- commandArgs(trailingOnly = TRUE)

f_rdata <- args[1]
out <- args[2]

load(f_rdata)

cust_func <- function() {
  MESuSiE::meSuSie_core(ld_list,ss_list,
    L=10, estimate_residual_variance =TRUE, max_iter =500)
}

benchmark_results1 <- microbenchmark(cust_func(), times=1, unit="second")

res <- tibble(method = "mesusie",
  trait = trait_name,
  "time" = summary(benchmark_results1)$mean)

write_tsv(res, out)
