library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

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
  
  xmap_res <- XMAP::XMAP_i(df_x1, df_x2, df_y1, df_y2, K = L2,
    Sigma1=c(0.001,0.001), Sigma2 = c(0.001,0.001), rho = 0.1,
    Omega = OmegaHat, Sig_E1 = 1, Sig_E2 = 1,
    prior_weights = rep(1/ncol(df_x1), ncol(df_x1)),
    maxIter = 500, tol = 1e-4,
    estimate_prior_variance = T,
    estimate_residual_variance = T,
    estimate_background_variance = F, initialize = F)
  
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

res1 <- NULL
for (jdx in 1:L2) {
  res1 <- cbind(res1, xmap_res$gamma[jdx] * xmap_res$mu1[jdx,])
}

res2 <- NULL
for (jdx in 1:L2) {
  res2 <- cbind(res2, xmap_res$gamma[jdx] * xmap_res$mu2[jdx,])
}

all_res <- cbind(rowSums(res1), rowSums(res2))
all_res <- tibble(as.data.frame(all_res))

write_tsv(all_res, out, col_names = FALSE)

