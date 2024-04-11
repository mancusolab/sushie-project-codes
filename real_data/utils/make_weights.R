library(tidyverse)
library(glue)

# proteins - sushie
study <- "mesa.proteins"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_proteins/weights",
  full.names = TRUE, pattern = key_name)
ct <- 1
for (lf in lfs) {
  print(ct)
  pref <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lf)
  df_tmp1 <- read_tsv(lf, col_types = cols())
  df_tmp2 <- read_tsv(glue("{pref}.normal.mega.weights.tsv"), col_types = cols())
  
  if (nrow(df_tmp1) != nrow(df_tmp2)) {
    stop()
  }
  
  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight")) %>%
    left_join(df_tmp2 %>%
        select(snp, mega_weight),
      by = "snp") %>%
    mutate(a11 = a1) %>%
    unite(col="rsid", chrom, pos, a1, a0, sep = ":") %>%
    mutate(rsid = paste0("chr", rsid)) %>%
    select(rsid, a11, contains("weight"))%>%
    rename(a1 = a11)
  
  colnames(new_df)[3:5] <- c("EUR_weight", "AFR_weight", "HIS_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  
  df_snps <- new_df %>%
    select(rsid)
  
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/all_of_us/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  
  write_tsv(df_snps,
    glue("/project/nmancuso_8/data/sushie/all_of_us/snps/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.snps.tsv"))
  ct <- ct + 1
}


# rnaseq - sushie
study <- "mesa.mrna"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_rnaseq/weights",
  full.names = TRUE, pattern = key_name)
ct <- 1
for (lf in lfs) {
  print(ct)
  pref <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lf)
  df_tmp1 <- read_tsv(lf, col_types = cols())
  df_tmp2 <- read_tsv(glue("{pref}.normal.mega.weights.tsv"),
    col_types = cols())
  
  if (nrow(df_tmp1) != nrow(df_tmp2)) {
    stop()
  }
  
  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight")) %>%
    left_join(df_tmp2 %>%
        select(snp, mega_weight),
      by = "snp") %>%
    mutate(a11 = a1) %>%
    unite(col="rsid", chrom, pos, a1, a0, sep = ":") %>%
    mutate(rsid = paste0("chr", rsid)) %>%
    select(rsid, a11, contains("weight")) %>%
    rename(a1 = a11)
  
  colnames(new_df)[3:5] <- c("EUR_weight", "AFR_weight", "HIS_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  
  df_snps <- new_df %>%
    select(rsid)
  
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/all_of_us/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  
  write_tsv(df_snps,
    glue("/project/nmancuso_8/data/sushie/all_of_us/snps/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.snps.tsv"))
  ct <- ct + 1
}


# genoa - sushie
study <- "genoa.mrna"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_genoa/weights",
  full.names = TRUE, pattern = key_name)
ct <- 1
for (lf in lfs) {
  print(ct)
  
  pref <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lf)
  df_tmp1 <- read_tsv(lf, col_types = cols())
  df_tmp2 <- read_tsv(glue("{pref}.normal.mega.weights.tsv"),
    col_types = cols())
  
  if (nrow(df_tmp1) != nrow(df_tmp2)) {
    stop()
  }
  
  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight")) %>%
    left_join(df_tmp2 %>%
        select(snp, mega_weight),
      by = "snp") %>%
    mutate(a11 = a1) %>%
    unite(col="rsid", chrom, pos, a1, a0, sep = ":") %>%
    mutate(rsid = paste0("chr", rsid)) %>%
    select(rsid, a11, contains("weight")) %>%
    rename(a1 = a11)
  
  colnames(new_df)[3:5] <- c("EUR_weight", "AFR_weight", "mega_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  
  df_snps <- new_df %>%
    select(rsid)
  
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/all_of_us/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  
  write_tsv(df_snps,
    glue("/project/nmancuso_8/data/sushie/all_of_us/snps/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.snps.tsv"))
  ct <- ct + 1
}

# 
# lfs <- list.files("/scratch1/zeyunlu/sushie_rnaseq/weights",
#   full.names = TRUE, pattern = key_name)
# 
# df_tmp <- tibble(full = lfs) %>%
#   mutate(gene = gsub("/scratch1/zeyunlu/sushie_rnaseq/weights/", "", gsub("\\.normal\\.sushie\\.weights\\.tsv", "", full)))
# 
# haha <- list.files("/project/nmancuso_8/data/sushie/weights/mesa.mrna",
#   recursive=TRUE, full.names = FALSE)
# 
# df_tmp2 <- tibble(full = haha) %>%
#   mutate(gene = gsub("chr[0-9]+/mesa.mrna.chr[0-9]+.", "", gsub("\\.sushie\\.weights\\.tsv", "", full)))
# 
# new_haha <- df_tmp %>%
#   filter(!gene %in% df_tmp2$gene)
# lfs <- new_haha$full


# make annot files
# proteins
lfs <- list.files("/scratch1/zeyunlu/sushie_proteins/weights",
  pattern = "normal.sushie")

lfs <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lfs)
df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_proteins_gene_list_noMHC.tsv",
  col_names = FALSE)

df_ref <- df_ref %>%
  filter(X15 %in% lfs)

write_tsv(df_ref, "/project/nmancuso_8/data/sushie/weights/mesa.proteins.genes.metadata.tsv", col_names = FALSE)

# rnaseq
lfs <- list.files("/scratch1/zeyunlu/sushie_rnaseq/weights",
  pattern = "normal.sushie")

lfs <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lfs)
df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_rnaseq_gene_list_noMHC.tsv",
  col_names = FALSE)

df_ref <- df_ref %>%
  filter(X12 %in% lfs)

write_tsv(df_ref, "/project/nmancuso_8/data/sushie/weights/mesa.mrna.genes.metadata.tsv", col_names = FALSE)


lfs <- list.files("/scratch1/zeyunlu/sushie_genoa/weights",
  pattern = "normal.sushie")

lfs <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lfs)
df_ref <- read_tsv("/scratch1/zeyunlu/sushie/genoa_sushie_gene_list_noMHC.tsv",
  col_names = FALSE)

df_ref <- df_ref %>%
  filter(X2 %in% lfs)

write_tsv(df_ref, "/project/nmancuso_8/data/sushie/weights/genoa.mrna.genes.metadata.tsv", col_names = FALSE)



# make credible set anno files
# proteins
proteins_lfs <- read_tsv("~/data/sushie/real/proteins_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait)

proteins_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_proteins_gene_list_noMHC.tsv",
  col_names = FALSE)

rnaseq_lfs <- read_tsv("~/data/sushie/real/rnaseq_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait)

rnaseq_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_rnaseq_gene_list_noMHC.tsv",
  col_names = FALSE)

genoa_lfs <- read_tsv("~/data/sushie/real/genoa_normal.sushie_cs.tsv.gz") %>%
  filter(!is.na(snp)) %>%
  distinct(trait)

genoa_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/genoa_sushie_gene_list_noMHC.tsv",
  col_names = FALSE)

res <- proteins_df_ref %>%
  filter(X15 %in% proteins_lfs$trait) %>%
  select(X1 = X4, X2 = X15) %>%
  mutate(X3 = "mesa.proteins") %>%
  bind_rows(rnaseq_df_ref %>%
      filter(X12 %in% rnaseq_lfs$trait) %>%
      select(X1 = X4, X2 = X12) %>%
      mutate(X3 = "mesa.mrna"),
    genoa_df_ref %>%
      filter(X2 %in% genoa_lfs$trait) %>%
      select(X1 = X4, X2) %>%
      mutate(X3 = "genoa.mrna"))

write_tsv(res, "/project/nmancuso_8/data/sushie/aou_csonly/metadata.tsv", col_names = FALSE)

for (idx in 1:22) {
  df_tmp <- res %>%
    filter(X1 == idx)
  write_tsv(df_tmp, glue("/project/nmancuso_8/data/sushie/aou_csonly/metadata_chr/metadata.chr{idx}.tsv"), col_names = FALSE)
}


traits <- c("BAS", "EOS", "HCT", "HGB", "LYM",
  "MCH", "MCHC", "MCV", "MON", "MPV", "NEU",
  "PLT", "RBC", "RDW", "WBC")

ans1 <- c("EUR", "AFR", "HIS")
ans2 <- c("EUR", "AFR")

study1 <- c("mesa.mrna", "mesa.proteins")
study2 <- c("genoa.mrna")

meta1 <- crossing(traits, ans1, study1) %>%
  rename(ans = ans1,
    study = study1) %>%
  bind_rows(crossing(traits, ans2, study2) %>%
      rename(ans = ans2,
        study = study2)) %>%
  left_join(res, by = c("study" = "X3"))

write_tsv(meta1, "/project/nmancuso_8/data/sushie/twas_metadata.tsv", col_names = FALSE)

test_meta <- meta1 %>%
  filter(study == "mesa.mrna") %>%
  filter(traits %in% c("RBC", "WBC")) %>%
  filter(X1 %in% 22) %>%
  filter(ans == "EUR")

test_meta <- test_meta[sample(nrow(test_meta), 40),]

write_tsv(test_meta, "/project/nmancuso_8/data/sushie/test_metadata.tsv", col_names = FALSE)


ans1 <- c("EUR", "AFR", "HIS")
ans2 <- c("EUR", "AFR")

study1 <- c("mesa.mrna", "mesa.proteins")
study2 <- c("genoa.mrna")

meta1 <- crossing(ans1, study1) %>%
  rename(ans = ans1,
    study = study1) %>%
  bind_rows(crossing( ans2, study2) %>%
      rename(ans = ans2,
        study = study2)) %>%
  left_join(res, by = c("study" = "X3"))

write_tsv(meta1, "/project/nmancuso_8/data/sushie/twas/new_metadata.tsv", col_names = FALSE)

for (idx in 1:22) {
  df_tmp <- meta1 %>%
    filter(X1 == idx)
  write_tsv(df_tmp, glue("/project/nmancuso_8/data/sushie/twas/metadata_chr/metadata.chr{idx}.tsv"), col_names = FALSE)
}

test_meta <- meta1 %>%
  filter(study == "mesa.mrna") %>%
  filter(X1 %in% 22) %>%
  filter(ans == "EUR")

test_meta <- test_meta[sample(nrow(test_meta), 100),]

write_tsv(test_meta, "/project/nmancuso_8/data/sushie/twas/new_test_metadata.tsv", col_names = FALSE)




###  add more
# make credible set anno files
proteins_lfs <- read_tsv("~/data/sushie/real/proteins_normal.sushie_cs.tsv.gz") %>%
  filter(is.na(snp)) %>%
  distinct(trait)

proteins_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_proteins_gene_list_noMHC.tsv",
  col_names = FALSE)

rnaseq_lfs <- read_tsv("~/data/sushie/real/rnaseq_normal.sushie_cs.tsv.gz") %>%
  filter(is.na(snp)) %>%
  distinct(trait)

rnaseq_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa_rnaseq_gene_list_noMHC.tsv",
  col_names = FALSE)

genoa_lfs <- read_tsv("~/data/sushie/real/genoa_normal.sushie_cs.tsv.gz") %>%
  filter(is.na(snp)) %>%
  distinct(trait)

genoa_df_ref <- read_tsv("/scratch1/zeyunlu/sushie/genoa_sushie_gene_list_noMHC.tsv",
  col_names = FALSE)

res <- proteins_df_ref %>%
  filter(X15 %in% proteins_lfs$trait) %>%
  select(X1 = X4, X2 = X15) %>%
  mutate(X3 = "mesa.proteins") %>%
  bind_rows(rnaseq_df_ref %>%
      filter(X12 %in% rnaseq_lfs$trait) %>%
      select(X1 = X4, X2 = X12) %>%
      mutate(X3 = "mesa.mrna"),
    genoa_df_ref %>%
      filter(X2 %in% genoa_lfs$trait) %>%
      select(X1 = X4, X2) %>%
      mutate(X3 = "genoa.mrna"))

write_tsv(res, "/project/nmancuso_8/data/sushie/aou_csonly/metadata_leftover.tsv", col_names = FALSE)

for (idx in 1:22) {
  df_tmp <- res %>%
    filter(X1 == idx)
  write_tsv(df_tmp, glue("/project/nmancuso_8/data/sushie/aou_csonly/metadata_chr_leftover/metadata.chr{idx}_leftover.tsv"), col_names = FALSE)
}



# new twas metadata
proteins_lfs <- read_tsv("~/Documents/github/data/sushie_results/real/proteins_normal.sushie_cs.tsv.gz") %>%
  distinct(trait)

proteins_df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_proteins_gene_list_noMHC.tsv",
  col_names = FALSE)

rnaseq_lfs <- read_tsv("~/Documents/github/data/sushie_results/real/rnaseq_normal.sushie_cs.tsv.gz") %>%
  distinct(trait)

rnaseq_df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_rnaseq_gene_list_noMHC.tsv",
  col_names = FALSE)

genoa_lfs <- read_tsv("~/Documents/github/data/sushie_results/real/genoa_normal.sushie_cs.tsv.gz") %>%
  distinct(trait)

genoa_df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/genoa_sushie_gene_list_noMHC.tsv",
  col_names = FALSE)

res <- proteins_df_ref %>%
  filter(X15 %in% proteins_lfs$trait) %>%
  select(X1 = X4, X2 = X15) %>%
  mutate(X3 = "mesa.proteins") %>%
  bind_rows(rnaseq_df_ref %>%
      filter(X12 %in% rnaseq_lfs$trait) %>%
      select(X1 = X4, X2 = X12) %>%
      mutate(X3 = "mesa.mrna"),
    genoa_df_ref %>%
      filter(X2 %in% genoa_lfs$trait) %>%
      select(X1 = X4, X2) %>%
      mutate(X3 = "genoa.mrna"))

for (idx in 1:22) {
  df_tmp <- res %>%
    filter(X1 == idx)
  write_tsv(df_tmp, glue("~/Downloads/twas_meta/metadata.chr{idx}.tsv"), col_names = FALSE)
}


df <- res %>%
  filter(X1 == 22) %>%
  group_by(X3) %>%
  slice_head(n = 1)

write_tsv(df, "~/Downloads/metadata.chr22.test.tsv", col_names = FALSE)

res %>%
  group_by(X1) %>%
  summarize(n = n())


df_fish <- read_tsv("~/Downloads/fish.tsv")
any(duplicated(df_fish))
all(df_fish$name %in% res$X2)

res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n()))) %>%
  group_by(X1, group) %>%
  distinct(X1, group)

res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n())),
    group = group + (X1 - 1)*22) %>%
  group_by(group) %>%
  distinct(group)

res2 <- res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n())),
    group = group + (X1-1)*22) %>%
  arrange(X1, group, X2)

for (idx in unique(res2$group)) {
  df_tmp <- res2 %>%
    filter(group == idx)
  write_tsv(df_tmp, glue("~/Downloads/twas_meta2/metadata.group{idx}.tsv"), col_names = FALSE)
}

all_res <- res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n())),
    group = group + (X1 - 1)*22) %>%
  group_by(X1, group) %>%
  distinct(X1, group) %>%
  arrange(X1, group)

write_tsv(all_res, "~/Downloads/meta_starter.tsv", col_names = FALSE)




# make v5 weights

# rnaseq - sushie
study <- "mesa.v5"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_v5/weights",
  full.names = TRUE, pattern = key_name)
# 2000
ct <- 20001
for (lf in lfs[20001:21695]) {
  print(ct)
  pref <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lf)
  df_tmp1 <- read_tsv(lf, col_types = cols())
  df_tmp2 <- read_tsv(glue("{pref}.normal.mega.weights.tsv"),
    col_types = cols())
  
  if (nrow(df_tmp1) != nrow(df_tmp2)) {
    stop()
  }
  
  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight")) %>%
    left_join(df_tmp2 %>%
        select(snp, mega_weight),
      by = "snp") %>%
    mutate(a11 = a1) %>%
    unite(col="rsid", chrom, pos, a1, a0, sep = ":") %>%
    mutate(rsid = paste0("chr", rsid)) %>%
    select(rsid, a11, contains("weight")) %>%
    rename(a1 = a11)
  
  colnames(new_df)[3:5] <- c("EUR_weight", "AFR_weight", "HIS_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  
  df_snps <- new_df %>%
    select(rsid)
  
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/aou_v5/weights/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  
  write_tsv(df_snps,
    glue("/project/nmancuso_8/data/sushie/aou_v5/snps/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.snps.tsv"))
  ct <- ct + 1
}

df_ref <- read_tsv("/scratch1/zeyunlu/sushie/mesa__gene_list_noMHC.tsv",
  col_names = FALSE)

lfs <- list.files("/scratch1/zeyunlu/sushie_v5/weights",
  pattern = key_name)

lfs_names <- gsub("\\.normal\\.sushie\\.weights\\.tsv", "", lfs)

df_meta <- df_ref %>%
  filter(X12 %in% lfs_names)

for(idx in 1:22) {
  tmp_meta <- df_meta %>%
    filter(X4 == idx)
  write_tsv(tmp_meta, glue("/project/nmancuso_8/data/sushie/aou_v5/metadata/metadata.chr{idx}.tsv"), col_names = FALSE)
}


# new twas metadata
v5_lfs <- read_tsv("~/Documents/github/data/sushie_results/real/v5_normal.sushie_cs.tsv.gz") %>%
  distinct(trait)

v5_df_ref <- read_tsv("~/Documents/github/data/sushie_results/metadata/mesa_rnaseq_v5_gene_list_noMHC.tsv",
  col_names = FALSE)

res <- v5_df_ref %>%
  filter(X12 %in% v5_lfs$trait) %>%
  select(X1 = X4, X2 = X12)


df_fish <- read_tsv("~/Downloads/fish_4.tsv")
any(duplicated(df_fish))
all(df_fish$name %in% res$X2)

res2 <- res %>%
  filter(!X2 %in% df_fish$name) %>%
  group_by(X1) %>%
  mutate(group = as.integer(gl(n(), 10, n())),
    group = group + (X1-1)*22) %>%
  arrange(X1, group, X2)

for (idx in unique(res2$group)) {
  df_tmp <- res2 %>%
    filter(group == idx)
  write_tsv(df_tmp, glue("~/Downloads/twas_metadata2/metadata.group{idx}.tsv"), col_names = FALSE)
}

all_res <- res2 %>%
  distinct(X1, group) %>%
  arrange(X1, group)

write_tsv(all_res, "~/Downloads/meta_starter.tsv", col_names = FALSE)



all_res %>% filter(group == 51)


# final

# proteins - sushie
study <- "mesa.proteins"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_proteins/weights",
  full.names = TRUE, pattern = key_name)
ct <- 1
for (lf in lfs) {
  print(ct)
  df_tmp1 <- read_tsv(lf, col_types = cols())

  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight"))
  
  colnames(new_df)[7:9] <- c("EUR_weight", "AFR_weight", "HIS_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/manuscript_data/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  
  ct <- ct + 1
}


# rnaseq - sushie
study <- "mesa.mrna"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_rnaseq/weights",
  full.names = TRUE, pattern = key_name)
ct <- 20001
for (lf in lfs[20001:22000]) {
  print(ct)
  df_tmp1 <- read_tsv(lf, col_types = cols())

  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight"))
  
  colnames(new_df)[7:9] <- c("EUR_weight", "AFR_weight", "HIS_weight")
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]
  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/manuscript_data/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))
  ct <- ct + 1
}


# genoa - sushie
study <- "genoa.mrna"
key_name <- "normal.sushie"
lfs <- list.files("/scratch1/zeyunlu/sushie_genoa/weights",
  full.names = TRUE, pattern = key_name)
ct <- 11200
for (lf in lfs[11200:12000]) {
  print(ct)
  df_tmp1 <- read_tsv(lf, col_types = cols())
  
  new_df <- df_tmp1 %>%
    select(trait, chrom, snp, pos, a0, a1, contains("weight"))
  
  colnames(new_df)[7:8] <- c("EUR_weight", "AFR_weight")
  
  gene <- df_tmp1$trait[1]
  chrom <- df_tmp1$chrom[1]

  write_tsv(new_df,
    glue("/project/nmancuso_8/data/sushie/manuscript_data/weights/{study}/chr{chrom}/{study}.chr{chrom}.{gene}.sushie.weights.tsv"))

  ct <- ct + 1
}


