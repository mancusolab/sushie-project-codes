# this is to randomly select 500 genes for simulaiton

library(tidyverse)

dd <- read_tsv("~/USCHPC/trash/gencode.v26lift37.GRCh37.genes.only.tsv") %>%
  filter(type %in% "protein_coding") %>%
  mutate(chr = as.numeric(gsub("chr", "", chr))) %>%
  filter(!is.na(chr)) %>%
  select(symbol, chr, tss, tes) %>%
  arrange(symbol, chr)

set.seed(1234)
idx <- sample(nrow(dd), nrow(dd), replace = FALSE)
new_dd <- dd[idx, ]
write_tsv(new_dd, file = "shuffle_genes_list.tsv", escape = "none", col_names = FALSE)

