library(tidyverse)
library(glue)
library(readxl)

dd <- read_excel("~/Downloads/41586_2021_3552_MOESM6_ESM.xlsx", skip = 2)

# prepare for liftover

df_final <- dd %>%
  filter(!is.na(`chrom.`)) %>%
  pivot_longer(cols=c(-`chrom.`, -start, -end)) %>%
  filter(value != 0) %>%
  group_by(`chrom.`, start, end) %>%
  mutate(name = gsub("[^a-zA-Z\\s]+", "-", name),
    name = gsub("\\s+", "-", name)) %>%
  summarize(new_name = paste(name, collapse = ",")) %>%
  mutate(`chrom.` = paste0("chr",`chrom.`))

# for liftover
write_tsv(df_final, "~/Documents/github/data/chiou_grch37.bed", col_names = FALSE)

# for reference
write_tsv(df_final, "~/Documents/github/data/chiou_et_al_grch37.bed.gz", col_names = FALSE)


# read in grch38

df_chiou <- read_tsv("~/Documents/github/data/chieou_grch38", col_names = FALSE)
colnames(df_chiou) <- c("CHR", "START", "END", "ANNOTATIONS")

ref1 <- read_tsv("~/USCHPC/data/gencode.v34.gene.only.tsv.gz")

ref2 <- read_tsv("~/USCHPC/data/gencode.v34.gene.only.tsv.gz") %>%
  filter(TYPE == "protein_coding")

all_tss <- tibble()
for (idx in 1:22) {
  print(idx)
  tmp_comp <- df_chiou%>%
    filter(CHR == glue("chr{idx}"))
  
  tmp_ref1 <- ref1 %>%
    filter(CHR == idx) %>%
    mutate(newTSS = ifelse(STRAND == "+", TSS, TES)) %>%
    select(CHR, newTSS)
  
  tmp_ref2 <- ref2 %>%
    filter(CHR == idx) %>%
    mutate(newTSS = ifelse(STRAND == "+", TSS, TES)) %>%
    select(CHR, newTSS)
  
  ans1 <- NULL
  ans2 <- NULL
  
  for (jdx in 1:nrow(tmp_comp)) {
      tmp_s <- tmp_comp[[jdx,"START"]]
      tmp_e <- tmp_comp[[jdx,"END"]]
      ans1 <- c(ans1, any(tmp_s <= tmp_ref1$newTSS  & 
          tmp_ref1$newTSS <= tmp_e) * 1)
      ans2 <- c(ans2, any(tmp_s <= tmp_ref2$newTSS  &
          tmp_ref2$newTSS <= tmp_e) * 1)
  }
  
  tmp_comp$tss_all <- ans1
  tmp_comp$tss_protein <- ans2
  
  all_tss <- all_tss %>%
    bind_rows(tmp_comp)
}

new_df <- all_tss %>%
  mutate(ANNOTATIONS = ifelse(tss_all == 0 & tss_protein == 0,
    ANNOTATIONS,
    ifelse(tss_all == 1 & tss_protein == 0,
      paste0(ANNOTATIONS, ",tss_all"),
      ifelse(tss_all == 0 & tss_protein == 1,
        paste0(ANNOTATIONS, ",tss_protein"),
        ifelse(tss_all == 1 & tss_protein == 1,
          paste0(ANNOTATIONS, ",tss_all,tss_protein"), NA))))) %>%
  select(CHR, START, END, ANNOTATIONS)

write_tsv(new_df, "~/Documents/github/data/chiou_et_al_grch38.bed.gz", col_names = FALSE)


# Satpathy et al.

dd <- read_table("~/Documents/github/data/GSE129785_scATAC-PBMCs-Fresh.peaks.txt.gz") %>%
  separate(col=Feature, into = c("CHR", "START", "END"), sep = "_")

write_tsv(dd, "~/Documents/github/data/fresh_peaks.bed", col_names = FALSE)


dd <- read_table("~/Documents/github/data/GSE129785_scATAC-PBMCs-Frozen.peaks.txt.gz") %>%
  separate(col=Feature, into = c("CHR", "START", "END"), sep = "_")

write_tsv(dd, "~/Documents/github/data/frozen_peaks.bed", col_names = FALSE)


dd <- read_table("~/Documents/github/data/GSE129785_scATAC-PBMCs-FrozenSort.peaks.txt.gz") %>%
  separate(col=Feature, into = c("CHR", "START", "END"), sep = "_")

write_tsv(dd, "~/Documents/github/data/frozensort_peaks.bed", col_names = FALSE)


dd <- read_tsv("~/Documents/github/data/grch38_fresh_peaks.bed", col_names = FALSE) %>%
  mutate(ANNO = "snATAC-seq-fresh-peaks")

write_tsv(dd, "~/Documents/github/data/snATAC-seq-fresh-peaks.bed.gz", col_names = FALSE)

dd <- read_tsv("~/Documents/github/data/grch38_frozen_peaks.bed", col_names = FALSE) %>%
  mutate(ANNO = "snATAC-seq-frozen-peaks")

write_tsv(dd, "~/Documents/github/data/snATAC-seq-frozen-peaks.bed.gz", col_names = FALSE)

dd <- read_tsv("~/Documents/github/data/grch38_frozensort_peaks.bed", col_names = FALSE) %>%
  mutate(ANNO = "snATAC-seq-frozen-sort-peaks")

write_tsv(dd, "~/Documents/github/data/snATAC-seq-frozen-sort-peaks.bed.gz", col_names = FALSE)





