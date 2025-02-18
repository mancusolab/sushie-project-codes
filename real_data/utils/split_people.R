library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

path_pt <- args[1]
pre_seed <- as.numeric(args[2])
path_out <- args[3]
set.seed(pre_seed)
# path_pt <- "/project/nmancuso_8/data/sushie/meta_data/mesa_rnaseq_pt_EUR.tsv"
df_pt <- read_tsv(path_pt, col_names = FALSE)

df_pt1 <- df_pt %>%
  mutate(group = sample(1:5, n(), replace = TRUE))

for (idx in 1:5) {
  tmp_df1 <- df_pt1 %>%
    filter(group == idx) %>%
    select(X1)
  tmp_df2 <- df_pt1 %>%
    filter(group != idx) %>%
    select(X1)
  
  write_tsv(tmp_df1, glue("{path_out}.test.pt.group{idx}.tsv"), col_names = FALSE)
  write_tsv(tmp_df2, glue("{path_out}.train.pt.group{idx}.tsv"), col_names = FALSE)
}
