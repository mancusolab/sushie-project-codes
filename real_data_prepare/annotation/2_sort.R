library(tidyverse)
library(glue)

setwd("/project/nmancuso_8/data/LDSC/baseline_annotation/bed_hg38")
lfs <- list.files()

notwant <- c("hg38_Human_Enhancer_Villar_Species_Enhancer_Count.bed",
  "hg38_GTEx_FE_META_TISSUE_GE_MaxCPP.bed",
  "hg38_BLUEPRINT_FE_META_TISSUE_H3K4me1_MaxCPP.bed",
  "hg38_BLUEPRINT_FE_META_TISSUE_H3K27ac_MaxCPP.bed",
  "hg38_BLUEPRINT_FE_META_TISSUE_DNAMETH_MaxCPP.bed",
  "hg38_Backgrd_Selection_Stat.bed",
  "hg38_ASMC.bed",
  "hg38_alleleage.bed")

lfs <- lfs[!lfs %in% notwant]
out <- "/scratch1/zeyunlu/genoa_enrich/annotation/"
for (lf in lfs) {
  tmp <- read_tsv(lf, col_names = FALSE) %>%
    arrange(X1, X2, X3) %>%
    mutate(X1 = as.numeric(gsub("chr", "", X1))) %>%
    filter(!is.na(X1))
  
  write_tsv(tmp, paste0(out, lf), col_names = FALSE)
}

setwd("/scratch1/zeyunlu/new_annot/raw")

lfs <- read_tsv("~/trash/ldsc_list.tsv", col_names=FALSE)

for (lf in lfs$X1) {
  new_name <- gsub("\\.", "_", gsub("\\.bed", "", gsub("hg38_", "LDSC_", lf)))
  df_tmp <- read_tsv(lf, col_names = FALSE) %>%
    filter(X2 < X3) %>%
    mutate(anno = new_name)
  
  write_tsv(df_tmp, glue("../original/{new_name}.bed.gz"), col_names = FALSE)
}

