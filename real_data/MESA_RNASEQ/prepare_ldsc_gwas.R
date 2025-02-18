library(tidyverse)
library(broom)
library(glue)

# prodcue GWAS files
traits <- c("WBC", "EOS", "MON", "LYM", "NEU", "BAS")

setwd("/project/nmancuso_8/data/Chen_ME_blood_GWAS/processed/new_gwas/EA")
for (trait in traits) {
  print(trait)
  lfs <- list.files(".", pattern = trait)
  df_tmp <- lfs %>% map_df(read_tsv, col_types = cols())
  write_tsv(df_tmp, glue("/project/nmancuso_8/data/sushie/chen_traits/EA_{trait}.tsv")) 
}
