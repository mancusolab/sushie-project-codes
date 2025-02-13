library(tidyverse)
library(glue)

scratch_path <- "/project/nmancuso_8/data/sushie/meta_data"

# get gene name

dd1 <- read_tsv(
  "/project/nmancuso_8/sgopalan/scripts/FUSION_GENOA/genoa_fusion_params.tsv",
  col_names = FALSE)

dd2 <- read_tsv(
  "/project/nmancuso_8/sgopalan/scripts/FUSION_GENOA/genoa_fusion_params2.tsv",
  col_names = FALSE)

dd <- bind_rows(dd1, dd2) %>%
  group_by(X5) %>%
  filter(n() != 1) %>%
  filter("ea" %in% X1 & "aa" %in% X1) %>%
  ungroup() %>%
  select(NAME = X5) %>%
  distinct()

gencode <- read_tsv(
  "/project/nmancuso_8/data/GENCODE/gencode.v34.gene.only.tsv.gz") %>%
  distinct(NAME,.keep_all = TRUE)


pos_dd <- gencode %>%
    filter(NAME %in% dd$NAME) %>%
    arrange(ID2)

write_tsv(pos_dd,
  "/project/nmancuso_8/data/GENOA/processed/expression/genoa_sushie_gene_list.tsv", 
  col_names = TRUE)

write_tsv(pos_dd, glue("{scratch_path}/genoa_sushie_gene_list.tsv"),
  col_names = FALSE)

write_tsv(filter(pos_dd, MHC == 0),
  "/project/nmancuso_8/data/GENOA/processed/expression/genoa_sushie_gene_list_noMHC.tsv",
  col_names = TRUE)

write_tsv(filter(pos_dd, MHC == 0),
  glue("{scratch_path}/genoa_sushie_gene_list_noMHC.tsv"),
  col_names = FALSE)


# get covar
covar1 <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/covar/ea_covars.tsv") %>%
  mutate(sex = as.numeric(sex == 2)) %>%
  select(-FID)

covar2 <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/covar/aa_covars.tsv") %>%
  mutate(sex = as.numeric(sex == 2)) %>%
  select(-FID)

# write_tsv(covar1, "/project/nmancuso_8/data/GENOA/processed/covariates/ea_covars_sushie.tsv.gz", col_names = FALSE)
# write_tsv(covar2, "/project/nmancuso_8/data/GENOA/processed/covariates/aa_covars_sushie.tsv.gz", col_names = FALSE)

write_tsv(covar1, glue("{scratch_path}/ea_covars_sushie.tsv.gz"), col_names = FALSE)
write_tsv(covar2, glue("{scratch_path}/aa_covars_sushie.tsv.gz"), col_names = FALSE)


# get pheno
ea_raw <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/ea/ea_GRCh38_expr.bed.gz")
aa_raw <- read_tsv("/project/nmancuso_8/data/GENOA/FUSION/expr/aa/aa_GRCh38_expr.bed.gz")

ea_raw$sum <- rowSums(abs(ea_raw[5:ncol(ea_raw)]))

ea_tmp <- ea_raw %>%
  filter(gene %in% pos_dd$NAME) %>%
  select(-`# chr`, -start, -end) %>%
  group_by(gene) %>%
  filter(sum == max(sum)) %>%
  arrange(gene) %>%
  select(-sum)

ea_pheno <- t(ea_tmp[,-1]) %>%
  as.data.frame() %>%
  rownames_to_column("IID")

colnames(ea_pheno) <- c("IID", ea_tmp$gene)

aa_raw$sum <- rowSums(abs(aa_raw[5:ncol(aa_raw)]))

aa_tmp <- aa_raw %>%
  filter(gene %in% pos_dd$NAME) %>%
  select(-`# chr`, -start, -end) %>%
  group_by(gene) %>%
  filter(sum == max(sum)) %>%
  arrange(gene) %>%
  select(-sum)

aa_pheno <- t(aa_tmp[,-1]) %>%
  as.data.frame() %>%
  rownames_to_column("IID")

colnames(aa_pheno) <- c("IID", aa_tmp$gene)

# # 14034
# write_tsv(ea_pheno, "/project/nmancuso_8/data/GENOA/processed/expression/ea_pheno_sushie.tsv.gz")
# write_tsv(aa_pheno, "/project/nmancuso_8/data/GENOA/processed/expression/aa_pheno_sushie.tsv.gz")

write_tsv(ea_pheno, glue("{scratch_path}/ea_pheno_sushie.tsv.gz"))
write_tsv(aa_pheno, glue("{scratch_path}/aa_pheno_sushie.tsv.gz"))

