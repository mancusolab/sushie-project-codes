library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

f_rdata <- args[1]
out <- args[2]
load(f_rdata)

# run XMAP ss

xmap_res <- tryCatch({
  fit_step1 <- XMAP::estimate_gc(
    data.frame(Z = df_ss1$Z[keep_idx], N = df_ss1$N[keep_idx]),
    data.frame(Z = df_ss2$Z[keep_idx], N = df_ss2$N[keep_idx]),
    df_ldscore$X1[keep_idx], df_ldscore$X2[keep_idx], df_ldscore$X3[keep_idx],
    reg_w1 = ld_w1[keep_idx], reg_w2 = ld_w2[keep_idx],
    reg_wx = sqrt(ld_w1[keep_idx] * ld_w2[keep_idx]),
    constrain_intercept = F)
  
  # stage 2 of LDSC: fix intercepts and estimate slopes
  fit_step2 <- XMAP::estimate_gc(
    data.frame(Z = df_ss1$Z, N = df_ss1$N),
    data.frame(Z = df_ss2$Z, N = df_ss2$N),
    df_ldscore$X1, df_ldscore$X2, df_ldscore$X3,
    reg_w1 = ld_w1, reg_w2 = ld_w2,
    reg_wx = sqrt(ld_w1 * ld_w2),
    constrain_intercept = T, fit_step1$tau1$coefs[1], fit_step1$tau2$coefs[1], fit_step1$theta$coefs[1])
  
  OmegaHat <- diag(c(fit_step2$tau1$coefs[2], fit_step2$tau2$coefs[2])) 
  OmegaHat[1, 2] <- fit_step2$theta$coefs[2] # co-heritability
  
  OmegaHat[lower.tri(OmegaHat)] <- OmegaHat[upper.tri(OmegaHat)]
  
  c1 <- fit_step2$tau1$coefs[1] 
  c2 <- fit_step2$tau2$coefs[1] 
  
  ar_sigma <- array(rep(c(1e-3, 1e-4, 1e-4, 1e-3), L2), dim = c(2, 2, L2))
  
  xmap_res <- XMAP::XMAP(simplify2array(list(df_ld1, df_ld2)),
    cbind(df_ss1$Z, df_ss2$Z),
    n = num_N, K = L2, Sigma = ar_sigma,
    Omega = OmegaHat, Sig_E = c(1, 1),
    maxIter = 500, tol = 1e-4,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    estimate_background_variance = FALSE
  )
  cs1 <- XMAP::get_CS(xmap_res, Xcorr= total_df_ld, coverage = 0.95, min_abs_corr = 0.5)
  pip1 <- XMAP::get_pip(xmap_res$gamma)[df_causal$`SNPIndex_1based`]
  xmap_res
}, error = function(e) {
  print("xmap cannot converge")
  q()
})

cs1 <- XMAP::get_CS(xmap_res, Xcorr= total_df_ld, coverage = 0.95, min_abs_corr = 0.5)
pip1 <- XMAP::get_pip(xmap_res$gamma)
new_causal <- df_causal
new_causal$pip <- pip1[df_causal$`SNPIndex_1based`]
new_causal$cali <- as.numeric(df_causal$`SNPIndex_1based` %in% cs1$cs)
new_causal$sim <- sim
new_causal$locus <- locus
new_causal$method <- "xmap"
new_causal$ancestry <- c(rep("1", L2/2), rep("2", L2/2))
write_tsv(new_causal, glue("{out}.pip.tsv"))

df_fdr <- tibble("method" = "xmap")
df_fdr$fdr_cs <- length(cs1$cs)
ct <- 0
if (length(cs1$cs) != 0) {
  for (idx in 1:length(cs1$cs)) {
    if (any(cs1$cs[idx] %in% df_causal$SNPIndex_1based)) {
      ct <- ct +1
    }
  }
}

df_fdr$as_in <- ct
df_fdr$fdr_high <- sum(pip1 > 0.95)
df_fdr$as_high <- sum(pip1 > 0.95 & 1:length(pip1) %in% df_causal$`SNPIndex_1based`)
df_fdr$sim <- sim
df_fdr$locus <- locus

write_tsv(df_fdr, glue("{out}.fdr.tsv"))

df_rho <- tibble("locus" = locus, "sim" = sim, CSIndex = 1:L2)
df_rho$est_rho <- NA
df_rho$frob <- NA
for (idx in 1:L2) {
  tmp_matrix <- xmap_res$Sigma[,,idx]
  tmp_rho <- tmp_matrix[1,2] / (sqrt(tmp_matrix[1,1]) * sqrt(tmp_matrix[2,2]))
  tmp_frob <- norm(tmp_matrix, type = "F")
  df_rho$est_rho[idx] <- tmp_rho
  df_rho$frob[idx] <- tmp_frob
}

df_rho <- df_rho %>%
  arrange(desc(frob)) %>%
  mutate(CSIndex = row_number())

write_tsv(df_rho, glue("{out}.rho.tsv"))
