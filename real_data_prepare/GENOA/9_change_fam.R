library(tidyverse)

# ea_id <- read_tsv("/project/nmancuso_8/data/GENOA/sushie/ea_373_pt.id", col_names = FALSE)
# aa_id <- read_tsv("/project/nmancuso_8/data/GENOA/sushie/aa_441_pt.id", col_names = FALSE)


ea_key <- read_csv("/project/nmancuso_8/data/GENOA/processed/keys/EA_dbGap_GEO.csv")
aa_key <- read_csv("/project/nmancuso_8/data/GENOA/processed/keys/AA_dbGap_GEO.csv")

setwd("/project/nmancuso_8/data/GENOA/processed/genotype/plink/annotated_dbsnp155")

for (idx in 1:22) {
  fname <- paste0("ea_chr", idx, "old.fam")
  dd <- read_tsv(fname, col_names = FALSE) %>%
    left_join(ea_key %>%
        filter(!is.na(dbgap_PIN2)) %>%
        select(dbgap_PIN2, GEO_ID1),
      by = c("X2" = "dbgap_PIN2")) %>%
    mutate(X2 = ifelse(!is.na(GEO_ID1), GEO_ID1, X2),
      X1 = X2) %>%
    select(-GEO_ID1)
  
  write_tsv(dd,  paste0("ea_chr", idx, ".fam"), col_names = FALSE)
  
  fname <- paste0("aa_chr", idx, "old.fam")
  dd <- read_tsv(fname, col_names = FALSE) %>%
    left_join(aa_key %>%
        filter(!is.na(dbgap_PIN2)) %>%
        select(dbgap_PIN2, GEO_ID2),
      by = c("X2" = "dbgap_PIN2")) %>%
    mutate(X2 = ifelse(!is.na(GEO_ID2), GEO_ID2, X2),
      X1 = X2) %>%
    select(-GEO_ID2)
  
  write_tsv(dd,  paste0("aa_chr", idx, ".fam"), col_names = FALSE)
  
}



