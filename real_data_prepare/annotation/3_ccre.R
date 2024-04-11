library(tidyverse)
# pbmc

pbmc_anno <- c("TSS", "activated-CD-T",
  "adaptive-NK", "classical-monocyte",
  # "cytotoxic-CD-T", "cytotoxic-NK",
   "megakaryocyte", "memory-B", "memory-CD-T",
  "naive-B", "naive-T", "non-classical-monocyte",
  "regulatory-T") 

write_tsv(tibble(pbmc_anno),
  "~/USCHPC/trash/pbmc_diabetes_anno.tsv", col_names = FALSE)

tcell_anno <- c("activated-CD-T", "cytotoxic-CD-T", "memory-CD-T",
  "naive-T", "regulatory-T") 
write_tsv(tibble(tcell_anno),
  "~/USCHPC/trash/tcell_diabetes_anno.tsv", col_names = FALSE)

mono_anno <- c("classical-monocyte", "non-classical-monocyte") 
write_tsv(tibble(mono_anno),
  "~/USCHPC/trash/mono_diabetes_anno.tsv", col_names = FALSE)


df_chiou <- read_tsv("~/Documents/github/data/chieou_grch38", col_names = FALSE)

ano1 <- unique(unlist(strsplit(df_chiou$ANNOTATIONS, ",")))

write_tsv(tibble(haha=ano1[!ano1 %in% pbmc_anno]), "~/USCHPC/trash/anno_exclude.tsv", col_names = FALSE)

