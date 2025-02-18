library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_pred/sim3_locus165/other.in.sim3.locus165.rdata"
f_rdata <- args[1]
out <- args[2]
load(f_rdata)

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
  res <- tibble("weights1" = rep(0, nrow(df_ss2)), "weights2" = rep(0, nrow(df_ss2)))
  write_tsv(res, out, col_names = FALSE)
  q()
})


all_res <- NULL
# number of pop
for (idx in 1:2) {
  # L
  tmp_res <- NULL
  for (jdx in 1:L2) {
    tmp_res <- cbind(tmp_res, xmap_res$gamma[jdx] * xmap_res$mu[idx, jdx,])
  }
  all_res <- cbind(all_res, rowSums(tmp_res))
}

all_res <- tibble(as.data.frame(all_res))

write_tsv(all_res, out, col_names = FALSE)

