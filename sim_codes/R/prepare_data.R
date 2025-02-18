library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

ss1 <- args[1]
ss2 <- args[2]
ld1 <- args[3]
ld2 <- args[4]
xy1 <- args[5]
xy2 <- args[6]
ldsc <- args[7]
causal <- args[8]
N <- args[9]
L2 <- as.numeric(args[10])
sim <- args[11]
locus <- args[12]
out <- args[13]

df_causal <- read_tsv(causal)
num_N <- as.numeric(unlist(stringr::str_split(N, ":")))
df_ss1 <- read_tsv(ss1) %>%
  mutate(SNP = paste0("SNP", pos),
    Z = beta / se,
    N = num_N[1]) %>%
  select(snp, SNP, Beta=beta, Se = se, Z, N, POS=pos) %>%
  column_to_rownames(var = "snp") %>%
  as.data.frame()

df_ss2 <- read_tsv(ss2) %>%
  mutate(SNP = paste0("SNP", pos),
    Z = beta / se,
    N = num_N[2]) %>%
  select(snp, SNP, Beta=beta, Se = se, Z, N, POS=pos) %>%
  column_to_rownames(var = "snp") %>%
  as.data.frame()
ss_list <- list("pop1" = df_ss1, "pop2" = df_ss2)

df_ld1 <- read_tsv(ld1, col_names = FALSE) %>%
  as.matrix()
# it might lose precision when pandas output this
new_matrix <- df_ld1
diag(new_matrix) <- 1
new_matrix[upper.tri(new_matrix)] <- df_ld1[upper.tri(df_ld1)]
new_matrix[lower.tri(new_matrix)] <- t(df_ld1)[lower.tri(new_matrix)]
df_ld1 <- new_matrix
colnames(df_ld1) <- df_ss1$SNP
rownames(df_ld1) <- rownames(df_ss1)

df_ld2 <- read_tsv(ld2, col_names = FALSE) %>%
  as.matrix()
new_matrix <- df_ld2
diag(new_matrix) <- 1
new_matrix[upper.tri(new_matrix)] <- df_ld2[upper.tri(df_ld2)]
new_matrix[lower.tri(new_matrix)] <- t(df_ld2)[lower.tri(new_matrix)]
df_ld2 <- new_matrix
colnames(df_ld2) <- df_ss2$SNP
rownames(df_ld2) <- rownames(df_ss2)
ld_list <- list("pop1" = df_ld1, "pop2" = df_ld2)


# run XMAP ss
df_ldscore <- read_tsv(ldsc, col_names = FALSE)
# remove extreme value
keep_idx <- which(df_ss1$Z^2 < 30 & df_ss2$Z^2 < 30)

ld_w1 <- 1 / sapply(df_ldscore$X1, function(x) max(x, 1))
ld_w2 <- 1 / sapply(df_ldscore$X2, function(x) max(x, 1))

df_xy1 <- read_tsv(xy1, col_names = FALSE)
df_xy2 <- read_tsv(xy2, col_names = FALSE)

df_x1 <- as.matrix(df_xy1[,2:ncol(df_xy1)])
df_x2 <- as.matrix(df_xy2[,2:ncol(df_xy2)])
df_y1 <- as.matrix(df_xy1[,1])
df_y2 <- as.matrix(df_xy2[,1])

total_df_ld <- num_N[1]/sum(num_N) * df_ld1 + num_N[2]/sum(num_N) * df_ld2

save(df_ss1, df_ss2, df_ld1, df_ld2, df_ldscore,
  df_x1, df_x2, df_y1, df_y2, ld_w1, ld_w2, ss_list, ld_list,
  num_N, L2, sim, locus, df_causal, keep_idx, total_df_ld,
  file = out)
