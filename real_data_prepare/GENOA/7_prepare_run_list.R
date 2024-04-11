library(tidyverse)

ans <- c("ea", "aa")
chr <- 1:22
write_tsv(crossing(ans, chr), "/project/nmancuso_8/data/GENOA/sushie/ans_chr_list.tsv",col_names = FALSE)

newchr <- paste0("chr", 1:22)
write_tsv(tibble(newchr, chr),
  "/project/nmancuso_8/data/GENOA/sushie/chr_convert.tsv", col_names = FALSE)

