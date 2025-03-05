# this script is to create main figure 3
library(tidyverse)
library(glue)
library(plotgardener)

# change the data folder to the zenodo-downloaded data folder
data_folder <- "~/Documents/github/data/sushie_results/real2"
metadata_folder <- "~/Documents/github/data/sushie_results/metadata2"

rnaseq_cov <- read_tsv(glue("{data_folder}/rnaseq_normal.sushie_cs.tsv.gz"))

rnaseq_meta <- read_tsv(glue("{data_folder}/rnaseq_normal.meta_cs.tsv.gz")) %>%
  filter(!is.na(snp)) %>%
  group_by(trait) %>%
  filter(n() > 3) %>%
  distinct(trait)

v5_cov <- read_tsv(glue("{data_folder}/v5_normal.sushie_cs.tsv.gz"))

rnaseq_susiex <- read_tsv(glue("{data_folder}/rnaseq_susiex_cs.tsv.gz")) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_mesusie <- read_tsv(glue("{data_folder}/rnaseq_mesusie_cs.tsv.gz")) %>%
  mutate(study = "mesa.mrna") %>%
  select(study, trait, CSIndex, snp, pip = pip_all)

rnaseq_xmap <- read_tsv(glue("{data_folder}/rnaseq_xmap_cs.tsv.gz"))

rnaseq_xmap %>%
  filter(grepl("URGCP", trait))

rnaseq_xmap %>%
  filter(grepl("SNHG5", trait))

df_tmp <- rnaseq_cov %>%
  filter(!is.na(snp)) %>%
  filter(!trait %in% rnaseq_susiex$trait) %>%
  filter(!trait %in% rnaseq_mesusie$trait) %>%
  group_by(trait) %>%
  filter(n() < 2) %>%
  filter(trait %in% rnaseq_meta$trait)

# annotation file downloaded from WASHU epigenome browser
# these are large files not included in zenodo
lfs <- list.files("~/Downloads/encode_anno/", full.names = TRUE)

sel_snp <- tibble()
for (idx in 1:nrow(df_tmp)) {
  tmp_snp_pos <- df_tmp$pos[idx]
  tmp_chr <- paste0("chr",df_tmp$chrom[idx])
  for (jdx in lfs) 
    tmp_chipseq <- readBigwig(jdx,
      chrom = tmp_chr, chromstart = tmp_snp_pos, chromend = tmp_snp_pos)
  if (nrow(tmp_chipseq) != 0){
    sel_snp <- sel_snp %>%
      bind_rows(tmp_chipseq %>%
          mutate(snp = df_tmp$snp[idx],
            pos = df_tmp$pos[idx],
            trait = df_tmp$trait[idx],
            chipseq = jdx))
  }
}

sel_gene <- sel_snp %>%
  filter(score >= 20) %>%
  distinct(trait) %>%
  filter(!grepl("URGCP", trait))

sel_snp2 <- sel_snp %>%
  filter(score >= 20) 

df_ref <- read_tsv(glue("{metadata_folder}/mesa_rnaseq_gene_list_noMHC.tsv"), col_names = FALSE) %>%
  filter(X12 %in% sel_gene$trait)

# write_tsv(df_ref, "~/Downloads/case_gene.tsv", col_names = FALSE)

lfs2 <- list.files("~/Downloads/eqtl/", full.names = TRUE)

df_z <- tibble()
df_z_all <- tibble()
for (idx in lfs2) {
  if (grepl("EUR", idx)) {
    pop_name <- "EUR"
  } else if (grepl("AFR", idx)) {
    pop_name <- "AFR"
  } else {
    pop_name <- "HIS"
  }
  df_z <- df_z %>%
    bind_rows(
      read_tsv(idx) %>%
        filter(ID %in% sel_snp2$snp) %>%
        select(snp = ID, Z = T_STAT) %>%
        mutate(pop = pop_name)
    )
  match <- sub(".*/(ENSG[0-9]+_\\w+)\\..*", "\\1", idx)
  df_z_all <- df_z_all %>%
    bind_rows(
      read_tsv(idx) %>%
        select(snp = ID, Z = T_STAT) %>%
        mutate(pop = pop_name,
          trait = match)
    )
}

df_total <- df_z %>%
  left_join(sel_snp2) %>%
  select(seqnames, snp, Z, pop, trait, score, chipseq) %>%
  group_by(trait) %>%
  filter(sum(abs(Z) > 1.96) == 3)


lfs <- list.files("~/Downloads/encode_anno2/", full.names = TRUE)

sel_snp3 <- tibble()
tmp_chr <- "chr6"
tmp_snp_pos <- 85678170
for (jdx in lfs) {
  sel_snp3 <- sel_snp3 %>%
    bind_rows(
      readBigwig(jdx,
        chrom = tmp_chr, chromstart = tmp_snp_pos, chromend = tmp_snp_pos) %>%
        mutate(chipseq = jdx)
    )
}

