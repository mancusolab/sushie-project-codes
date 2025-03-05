library(tidyverse)
library(glue)
args <- commandArgs(trailingOnly = TRUE)

# prefix <- "/scratch1/zeyunlu/tmp_2pop/sim3_locus233"
# sim <- 3
# locus <- 233
# causal_file <- "/scratch1/zeyunlu/tmp_2pop/sim3_locus233/other.causal.sim3.locus233.tsv"
# L2 <- 2
# out <- "/scratch1/zeyunlu/sushie_sim_2pop/"
prefix <- args[1]
sim <- args[2]
locus <- args[3]
causal_file <- args[4]
out <- args[5]
L2 <- as.numeric(args[6])


# rho
# sushie
df_rho_sushie <- read_tsv(glue("{prefix}/other.sushie.sim{sim}.locus{locus}.rho.tsv"))

# xmap
xmap_rho_path <- glue("{prefix}/other.xmap.in.sim{sim}.locus{locus}.rho.tsv")

if (file.exists(xmap_rho_path)) {
  xmap_rho <- read_tsv(xmap_rho_path)
  df_rho_sushie <- df_rho_sushie %>%
    bind_rows(
      df_rho_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        left_join(xmap_rho %>%
            rename(Lidx = CSIndex) %>%
            mutate(method = "xmap") %>%
            select(-frob)
        )
    )
}

# xmap_ind
xmapind_rho_path <- glue("{prefix}/other.xmap.ind.sim{sim}.locus{locus}.rho.tsv")

if (file.exists(xmapind_rho_path)) {
  xmapind_rho <- read_tsv(xmapind_rho_path)
  df_rho_sushie <- df_rho_sushie %>%
    bind_rows(
      df_rho_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        left_join(xmapind_rho %>%
            rename(Lidx = CSIndex) %>%
            mutate(method = "xmap_ind") %>%
            select(-frob)
        )
    )
}

# mesusie
mesusie_rho_path <- glue("{prefix}/other.mesusie.in.sim{sim}.locus{locus}.rho.tsv")

if (file.exists(mesusie_rho_path)) {
  mesusie_rho <- read_tsv(mesusie_rho_path)
  df_rho_sushie <- df_rho_sushie %>%
    bind_rows(
      df_rho_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        left_join(mesusie_rho %>%
            rename(Lidx = CSIndex) %>%
            mutate(method = "mesusie") %>%
            select(-frob)
        )
    )
}

write_tsv(df_rho_sushie, glue("{out}/sim{sim}.locus{locus}.rho.tsv"))

# # susiex
# susiex_cs <- glue("{prefix}/susiex.in.sim{sim}.locus{locus}.cs")
# susiex_snp <- glue("{prefix}/susiex.in.sim{sim}.locus{locus}.snp")

# pip
# sushie
df_pip_sushie <- read_tsv(glue("{prefix}/other.sushie.sim{sim}.locus{locus}.pip.tsv"))

# susiex
if (file.exists(susiex_cs) & file.exists(susiex_snp)) {
  df_susiex_snp <- read_tsv(susiex_snp)
  n_cs <- sum(grepl("PIP", colnames(df_susiex_snp)))
  if (n_cs < L2) {
    for (idx in 1:(L2-n_cs)) {
      df_susiex_snp[[glue("PIP{idx}")]] <- rep(1/nrow(df_susiex_snp),
        nrow(df_susiex_snp))
    }
  }
  df_mat <- df_susiex_snp %>%
    select(contains("PIP")) %>%
    mutate(across(everything(), ~ 1 - .)) %>%
    rowwise() %>%
    mutate(result = 1 - prod(c_across(contains("PIP")))) %>% 
    ungroup()
  
  df_susiex_snp$susiex_pip <- df_mat$result 
  
  df_pip_sushie <- df_pip_sushie %>%
    left_join(df_susiex_snp %>%
        select(SNP, susiex_pip), by = c("snp" = "SNP"))
  
} else {
  df_pip_sushie$susiex_pip <- NA
}

# mesusie
mesusie_pip_path <- glue("{prefix}/other.mesusie.in.sim{sim}.locus{locus}.pip.tsv")
if (file.exists(mesusie_pip_path)) {
  df_mesusie_pip <- read_tsv(mesusie_pip_path)
  df_pip_sushie <- df_pip_sushie %>%
    left_join(df_mesusie_pip %>%
        rename(mesusie_pip = mesusie))
} else {
  df_pip_sushie$mesusie_pip <- NA
}

# xmap
xmap_pip_path <- glue("{prefix}/other.xmap.in.sim{sim}.locus{locus}.pip.tsv")
if (file.exists(xmap_pip_path)) {
  df_xmap_pip <- read_tsv(xmap_pip_path)
  df_pip_sushie <- df_pip_sushie %>%
    left_join(df_xmap_pip %>%
        rename(xmap_pip = xmap))
} else {
  df_pip_sushie$xmap_pip <- NA
}

# xmap ind
xmapind_pip_path <- glue("{prefix}/other.xmap.ind.sim{sim}.locus{locus}.pip.tsv")
if (file.exists(xmapind_pip_path)) {
  df_xmapind_pip <- read_tsv(xmapind_pip_path)
  df_pip_sushie <- df_pip_sushie %>%
    left_join(df_xmapind_pip %>%
        rename(xmap_ind_pip = xmap))
} else {
  df_pip_sushie$xmap_ind_pip <- NA
}

# cs
# sushie
df_cs_sushie <- read_tsv(glue("{prefix}/other.sushie.sim{sim}.locus{locus}.cs.tsv"))

# susiex
df_ref <- read_tsv(causal_file)
if (file.exists(susiex_cs)) {
  df_susiex_cs <- read_tsv(susiex_cs, col_names = FALSE)
  if (!df_susiex_cs$X1[1] %in% c("NULL", "FAIL")) {
    df_susiex_cs <- read_tsv(susiex_cs)
    df_cs_sushie <- df_cs_sushie %>%
      bind_rows(
        df_cs_sushie %>%
          distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
          mutate(method = "susiex") %>%
          left_join(df_susiex_cs %>%
              mutate(CSIndex = CS_ID,
                snp = SNP,
                method = "susiex",
                causal = ifelse(snp %in% df_ref$SNP, 1, 0)) %>%
              left_join(df_pip_sushie %>%
                  select(snp, pip_all = susiex_pip), by = "snp") %>%
              select(snp, CSIndex, causal, method, pip_all))
      )
  } else {
    df_cs_sushie <- df_cs_sushie %>%
      bind_rows(
        df_cs_sushie %>%
          distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
          mutate(method = "susiex") %>%
          left_join(tibble(
            method = "susiex"))
      )
  }
} else {
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "susiex") %>%
        left_join(tibble(
          method = "susiex"))
    )
}

# mesusie
mesusie_cs_path <- glue("{prefix}/other.mesusie.in.sim{sim}.locus{locus}.cs.tsv")

if (file.exists(mesusie_cs_path)) {
  df_mesusie_cs <- read_tsv(mesusie_cs_path) %>%
    left_join(df_pip_sushie %>%
        select(snp, SNPIndex_1based), by = "SNPIndex_1based") %>%
    mutate(causal = ifelse(snp %in% df_ref$SNP, 1, 0)) %>%
    left_join(df_pip_sushie %>%
        select(SNPIndex_1based, pip_all = mesusie_pip), by = "SNPIndex_1based")
  
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "mesusie") %>%
        left_join(df_mesusie_cs %>%
            mutate(method = "mesusie") %>%
            select(snp, CSIndex, SNPIndex_1based, causal, method, pip_all)
        )
    )
} else {
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "mesusie") %>%
        left_join(tibble(
          method = "mesusie"))
    )
}

# xmap
xmap_cs_path <- glue("{prefix}/other.xmap.in.sim{sim}.locus{locus}.cs.tsv")

if (file.exists(xmap_cs_path)) {
  df_xmap_cs <- read_tsv(xmap_cs_path) %>%
    left_join(df_pip_sushie %>%
        select(snp, SNPIndex_1based), by = "SNPIndex_1based") %>%
    mutate(causal = ifelse(snp %in% df_ref$SNP, 1, 0)) %>%
    left_join(df_pip_sushie %>%
        select(SNPIndex_1based, pip_all = xmap_pip), by = "SNPIndex_1based")
  
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "xmap") %>%
        left_join(df_xmap_cs %>%
            mutate(method = "xmap") %>%
            select(snp, CSIndex, causal, method, pip_all)
        )
    )
} else {
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "xmap") %>%
        left_join(tibble(
          method = "xmap"))
    )
}

# xmap ind
xmapind_cs_path <- glue("{prefix}/other.xmap.ind.sim{sim}.locus{locus}.cs.tsv")

if (file.exists(xmapind_cs_path)) {
  df_xmapind_cs <- read_tsv(xmapind_cs_path)  %>%
    left_join(df_pip_sushie %>%
        select(snp, SNPIndex_1based), by = "SNPIndex_1based") %>%
    mutate(causal = ifelse(snp %in% df_ref$SNP, 1, 0)) %>%
    left_join(df_pip_sushie %>%
        select(SNPIndex_1based, pip_all = xmap_ind_pip), by = "SNPIndex_1based")
  
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "xmap_ind") %>%
        left_join(df_xmapind_cs %>%
            mutate(method = "xmap_ind") %>%
            select(snp, CSIndex, SNPIndex_1based, causal, method, pip_all)
        )
    )
} else {
  df_cs_sushie <- df_cs_sushie %>%
    bind_rows(
      df_cs_sushie %>%
        distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
        mutate(method = "xmap_ind") %>%
        left_join(tibble(
          method = "xmap_ind"))
    )
}

df_cs_sushie <- df_pip_sushie %>%
  select(snp, SNPIndex_1based) %>%
  right_join(df_cs_sushie %>%
  select(-SNPIndex_1based, -SNPIndex_0based), by = "snp") %>%
  select(-snp)

write_tsv(df_cs_sushie, glue("{out}/sim{sim}.locus{locus}.cs.tsv"))

df_pip_sushie <- df_pip_sushie %>%
  select(-snp, -SNPIndex_0based)

write_tsv(df_pip_sushie, glue("{out}/sim{sim}.locus{locus}.pip.tsv"))

