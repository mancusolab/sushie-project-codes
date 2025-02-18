library(tidyverse)
library(glue)


# make annot files
lfs <- list.files("/project/nmancuso_8/data/sushie/aou_reg/snps", recursive = TRUE)

total_gene <- tibble(filepath = lfs) %>%
  mutate(filepath = gsub("\\.snps\\.tsv", "", filepath)) %>%
  separate(filepath, into = c("folder", "rest"), sep = "/", remove = TRUE) %>%
  separate(rest, into = c("type", "rest"), sep = "\\.chr[0-9]+\\.", remove = TRUE) %>%
  mutate(folder = as.numeric(gsub("chr", "", folder))) %>%
  arrange(folder, rest)

for (idx in 1:22) {
  tmp_df <- total_gene %>%
    filter(folder == idx)
  write_tsv(tmp_df, glue("/project/nmancuso_8/data/sushie/aou_reg/metadata/all.chr{idx}.metadata.tsv"), col_names = FALSE)
}

new_total <- total_gene %>%
  mutate(group = as.integer(gl(n(), 300, n())))

for (idx in unique(new_total$group)) {
  tmp_new <- new_total %>%
    filter(group == idx)
  write_tsv(tmp_new, glue("/project/nmancuso_8/data/sushie/aou_reg/metadata_twas/twas.group{idx}.meta.tsv"), col_names = FALSE)
}

# get starter

all_res <- res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n())),
    group = group + (X1 - 1)*22) %>%
  group_by(X1, group) %>%
  distinct(X1, group) %>%
  arrange(X1, group)

write_tsv(all_res, "~/Downloads/meta_starter.tsv", col_names = FALSE)




