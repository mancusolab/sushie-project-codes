library(tidyverse)
library(broom)
library(glue)

# create annotions
# sushie
base_path <- "/scratch1/zeyunlu/sushie_rnaseq/"
setwd(glue("{base_path}sushie/weights"))

lfs <- list.files(".", pattern = "normal\\.sushie")

all_res <- tibble(chrom=numeric(), snp=character(),
  pos=numeric(), a0=character(), a1=character(),
  sushie_alpha=numeric(), sushie_cs=numeric(),
  indep_alpha=numeric(), indep_cs=numeric(),
  meta_alpha=numeric(), meta_cs=numeric(),
  mega_alpha=numeric(), mega_cs=numeric(),
  mesusie_alpha=numeric(), mesusie_cs=numeric(),
  susiex_alpha=numeric(), susiex_cs=numeric())

ct <- 0

for (lf in lfs) {
  ct <- ct + 1
  print(ct)
  df_ref <- read_tsv(lf, col_types = cols(a1 = col_character(),
    a0 = col_character())) %>%
    select(chrom, snp, pos, a0, a1)
  old_nrow <- nrow(df_ref)
  tmp_name <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lf)
  
  # sushie
  tmp_cs <- read_tsv(
    glue("{base_path}sushie/cs/{tmp_name}.normal.sushie.cs.tsv"),
    col_types = cols(a1 = col_character(),
      a0 = col_character()))
  
  if (!is.na(tmp_cs$snp[1])) {
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, sushie_alpha = alpha) %>%
          mutate(sushie_cs = 1) %>%
          group_by(snp) %>%
          mutate(sushie_alpha = max(sushie_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  # indep
  tmp_cs <- read_tsv(
    glue("{base_path}sushie/cs/{tmp_name}.indep.sushie.cs.tsv"),
    col_types = cols(a1 = col_character(),
      a0 = col_character()))
  
  if (!is.na(tmp_cs$snp[1])) {
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, indep_alpha = alpha) %>%
          mutate(indep_cs = 1) %>%
          group_by(snp) %>%
          mutate(indep_alpha = max(indep_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  # meta
  tmp_cs <- read_tsv(
    glue("{base_path}sushie/cs/{tmp_name}.normal.meta.cs.tsv"),
    col_types = cols(a1 = col_character(),
      a0 = col_character()))
  
  if (!is.na(tmp_cs$snp[1])) {
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, meta_alpha = alpha) %>%
          mutate(meta_cs = 1) %>%
          group_by(snp) %>%
          mutate(meta_alpha = max(meta_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  # mega
  tmp_cs <- read_tsv(
    glue("{base_path}sushie/cs/{tmp_name}.normal.mega.cs.tsv"),
    col_types = cols(a1 = col_character(),
      a0 = col_character()))
  
  if (!is.na(tmp_cs$snp[1])) {
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, mega_alpha = alpha) %>%
          mutate(mega_cs = 1) %>%
          group_by(snp) %>%
          mutate(mega_alpha = max(mega_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  # mesusie
  file_path <- glue("{base_path}mesusie/cs/mesusie.{tmp_name}.cs.tsv")
  if (file.exists(file_path)) {
    tmp_cs <- read_tsv(file_path,
      col_types = cols(a1 = col_character(),
        a0 = col_character()))
    
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, mesusie_alpha = alpha) %>%
          mutate(mesusie_cs = 1) %>%
          group_by(snp) %>%
          mutate(mesusie_alpha = max(mesusie_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  # susiex
  file_path <- glue("{base_path}susiex/cs/susiex.{tmp_name}.cs.tsv")
  if (file.exists(file_path)) {
    tmp_cs <- read_tsv(file_path, col_types = cols())
    
    df_ref <- df_ref %>%
      left_join(tmp_cs %>%
          select(snp, susiex_alpha = alpha) %>%
          mutate(susiex_cs = 1) %>%
          group_by(snp) %>%
          mutate(susiex_alpha = max(susiex_alpha)) %>%
          distinct(snp, .keep_all = TRUE),
        by = "snp"
      )
  }
  
  if (ncol(df_ref) > 5){
    df_ref[is.na(df_ref)] <- 0
    df_ref <- df_ref[rowSums(df_ref[,6:ncol(df_ref)]) != 0,]
    all_res <- bind_rows(all_res, df_ref)
  }
}
all_res[is.na(all_res)] <- 0

all_res2 <- all_res %>%
  group_by(chrom, snp, pos, a0, a1) %>%
  summarize(sushie_alpha = max(sushie_alpha),
    sushie_cs = max(sushie_cs),
    indep_alpha = max(indep_alpha),
    indep_cs = max(indep_cs),
    meta_alpha = max(meta_alpha),
    meta_cs = max(meta_cs),
    mega_alpha = max(mega_alpha),
    mega_cs = max(mega_cs),
    mesusie_alpha = max(mesusie_alpha),
    mesusie_cs = max(mesusie_cs),
    susiex_alpha = max(susiex_alpha),
    susiex_cs = max(susiex_cs))

write_tsv(all_res2, "/project/nmancuso_8/data/sushie/annotations/mesa.rnaseq.sushie.annot.tsv.gz")
