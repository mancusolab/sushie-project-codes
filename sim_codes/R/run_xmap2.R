library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

f_rdata <- args[1]
out <- args[2]
out2 <- args[3]
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
    Sigma1=c(0.001, 0.001), Sigma2 = c(0.001,0.001), rho = 0.1,
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
  df_sens <- tibble("xmap" = 0)
  df_sens$sim <- sim
  df_sens$locus <- locus
  write_tsv(df_sens, glue("{out}.sens.tsv"))
  q()
})


cs1 <- XMAP::get_CS(xmap_res, Xcorr= total_df_ld, coverage = 0.95, min_abs_corr = 0.5)
pip1 <- XMAP::get_pip(xmap_res$gamma)[df_causal$`SNPIndex_1based`]

df_sens <- tibble("xmap" = 1)
df_sens$sim <- sim
df_sens$locus <- locus
write_tsv(df_sens, glue("{out}.sens.tsv"))

df_causal$xmap <- pip1
df_causal$sim <- sim
df_causal$locus <- locus
df_causal$cali <- as.numeric(df_causal$`SNPIndex_1based` %in% unlist(cs1$cs))
write_tsv(df_causal, glue("{out}.pip.tsv"))


df_cs <- tibble("CSIndex" = seq(L2))
df_cs$sim <- sim
df_cs$locus <- locus
df_cs$xmap <- NA

if (length(cs1$cs) != 0 ) {
  for (jdx in 1:length(cs1$cs)) {
    df_cs$xmap[jdx] <- length(cs1$cs[[jdx]])  
  }
}

write_tsv(df_cs, glue("{out}.cs.tsv"))

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


# additional

pip1 <- XMAP::get_pip(xmap_res$gamma)
df_pip <- tibble(SNPIndex_1based = 1:length(pip1), xmap = pip1)
df_pip$sim <- sim
df_pip$locus <- locus
write_tsv(df_pip, glue("{out2}.pip.tsv"))

if (length(cs1$cs) != 0 ) {
  df_cs <- tibble()
  for (jdx in 1:length(cs1$cs)) {
    df_cs <- df_cs %>%
      bind_rows(
        tibble(CSIndex = jdx, 
          SNPIndex_1based = cs1$cs[[jdx]]
        )
      )
  }
  
  df_cs$sim <- sim
  df_cs$locus <- locus
  
  write_tsv(df_cs, glue("{out2}.cs.tsv"))
}

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

write_tsv(df_rho, glue("{out2}.rho.tsv"))

