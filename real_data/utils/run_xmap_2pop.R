library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)
# f_rdata <- "/scratch1/zeyunlu/tmp_2pop/sim4_locus232/other.in.sim4.locus232.rdata"
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
  
  if (any(is.nan(OmegaHat))) {
    stop("OmegaHat contains nan")
  }
  
  c1 <- fit_step2$tau1$coefs[1] 
  c2 <- fit_step2$tau2$coefs[1] 
  
  ar_sigma <- array(rep(c(1e-3, 1e-4, 1e-4, 1e-3), 10), dim = c(2, 2, 10))
  
  xmap_res <- XMAP::XMAP(simplify2array(list(df_ld1, df_ld2)),
    cbind(df_ss1$Z, df_ss2$Z),
    n = num_N, K = 10, Sigma = ar_sigma,
    Omega = OmegaHat, Sig_E = c(1, 1),
    maxIter = 500, tol = 1e-4,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    estimate_background_variance = FALSE
  )
  
  cs1 <- XMAP::get_CS(xmap_res, Xcorr= total_df_ld, coverage = 0.95, min_abs_corr = 0.5)
  pip1 <- XMAP::get_pip(xmap_res$gamma)
  
  xmap_res
  
}, error = function(e) {
  print("xmap cannot converge")
  q()
})

cs1 <- XMAP::get_CS(xmap_res, Xcorr= total_df_ld, coverage = 0.95, min_abs_corr = 0.5)
pip1 <- XMAP::get_pip(xmap_res$gamma)

all_res <- NULL
for (idx in 1:2) {
  # L
  tmp_res <- NULL
  for (jdx in 1:10) {
    tmp_res <- cbind(tmp_res, xmap_res$gamma[jdx] * xmap_res$mu[idx, jdx,])
  }
  all_res <- cbind(all_res, rowSums(tmp_res))
}

df_ref$pop1_weights <- all_res[,1]
df_ref$pop2_weights <- all_res[,2]
df_ref$pip_all <- pip1
df_ref$trait <- trait_name
df_ref$method <- "xmap"
write_tsv(df_ref, glue("{out}/weights/xmap.{trait_name}.weights.tsv"))


df_cs <- tibble()
if (length(cs1$cs) != 0 ) {
  for (jdx in 1:length(cs1$cs)) {
    df_tmp <- tibble(SNPIndex = cs1$cs[[jdx]], CSIndex = jdx)
    df_cs <- bind_rows(df_cs, df_tmp)
  }
}

if (nrow(df_cs) != 0) {
  df_res <- df_cs %>%
    left_join(df_ref %>%
        select(chrom, snp, pos, a0, a1, pip_all, trait) %>%
        rownames_to_column() %>%
        rename(SNPIndex = rowname) %>%
        mutate(SNPIndex = as.numeric(SNPIndex)),
      by = "SNPIndex") %>%
    select(chrom, snp, pos, a0, a1, CSIndex, SNPIndex, pip_all, trait) %>%
    mutate(method = "xmap")
  
  write_tsv(df_res, glue("{out}/cs/xmap.{trait_name}.cs.tsv"))
  
  
  df_rho <- tibble()
  for (idx in 1:10) {
    tmp_matrix <- xmap_res$Sigma[,,idx]
    tmp_rho <- tmp_matrix[1,2] / (sqrt(tmp_matrix[1,1]) * sqrt(tmp_matrix[2,2]))
    tmp_frob <- norm(tmp_matrix, type = "F")
    df_rho <- df_rho %>%
      bind_rows(
        tibble(type = "ancestry1_ancestry2_est_corr",
          est_corr = tmp_rho,
          frob = tmp_frob),
      )
  }
  
  df_rho <- df_rho %>%
    arrange(desc(frob)) %>%
    mutate(CSIndex = row_number()) %>%
    mutate(trait = trait_name,
      method = "xmap")
  
  write_tsv(df_rho, glue("{out}/corr/xmap.{trait_name}.corr.tsv"))
}