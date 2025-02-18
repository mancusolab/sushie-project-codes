library(tidyverse)
library(glue)

args <- commandArgs(TRUE)
# weight_path <- "/project/nmancuso_8/data/sushie/aou_reg/weights/chr5/mesa.mrna.chr5.ENSG00000081189_MEF2C.weights.tsv"
# heri_path <- "/scratch1/zeyunlu/sushie_rnaseq/sushie/her/ENSG00000081189_MEF2C.normal.sushie.her.tsv"
# cvr2_path <- "/scratch1/zeyunlu/temp_cv_rnaseq/tempf_1/ENSG00000081189_MEF2C/ENSG00000081189_MEF2C.cvr2.tsv"


weight_path <- args[1]
heri_path <- args[2]
cvr2_path <- args[3]
trait_name <- args[4]
N_total <- args[5]
out_path <- args[6]

weight <- read_tsv(weight_path)

heri <- read_tsv(heri_path)
cvr2 <- read_tsv(cvr2_path) %>%
  filter(group == "match") %>%
  select(sushie, susie, enet, lasso, gblup = ridge, type, pop)

pop_list <- c("EUR", "AFR", "HIS")
total_ss <- as.numeric(unlist(strsplit(N_total, ":")))
for (idx in 1:3) {
  tmp_cvr2 <- cvr2 %>%
    filter(pop == pop_list[idx])
  N.tot <- total_ss[idx]
  
  # cv.performance
  cv.performance <- as.matrix(tmp_cvr2[,1:5])
  rownames(cv.performance) <- c("rsq", "pval")
  colnames(cv.performance) <- c("sushie", "susie", "enet", "lasso", "gblup")

  # hsq
  hsq <- c(as.numeric(heri[idx, 3]), NA)

  # hsq.pv
  hsq.pv <- heri[idx,5]

  # snps
  weight$cm <- 0
  snps <- weight[c("chrom", "snp", "cm", "pos", "a0", "a1")]
  colnames(snps) <- paste0("V", 1:6)
  snps <- data.frame(snps)

  # wgt.matrix
  tmp_weight <- weight[c(glue("sushie_pop{idx}"), "mega_weight",
    "enet_weight", "lasso_weight", "gblup_weight")]
  wgt.matrix <- as.matrix(tmp_weight)
  rownames(wgt.matrix) <- weight$snp
  colnames(wgt.matrix) <- c("sushie", "susie", "enet", "lasso", "gblup")
  
  # save r data
  save(cv.performance, hsq, hsq.pv, N.tot, snps, wgt.matrix,
    file = glue("{out_path}/{pop_list[idx]}/rdata/{trait_name}.{pop_list[idx]}.fusion.RData"))
  write_tsv(tibble(cbind(trait_name, pop_list[idx], hsq[1], hsq.pv)),
    glue("{out_path}/{pop_list[idx]}/hsq/{trait_name}.{pop_list[idx]}.hsq"),
    col_names = FALSE)
}
