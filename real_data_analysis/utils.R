library(tidyverse)
library(broom)
library(glue)

# color
# rnaseq, genoa, proteins
main_study_color <- c("#8dd3c7", "#fb8072", "#bebada")

ffont <- "sans"
fontsize <- 9
errbar_width <- 0.5
scaleFUN <- function(x) sprintf("%.2f", x)

p_width <- 8
p_height <- 2.5
point_size <- 1.5

theme_qtl <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    text=element_text(size = 14))
}


theme_tss <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    text=element_text(size = 14))
}

theme_valid <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    text=element_text(size = 14))
}


basic_sum <- function(..., total, subsets=NULL) {
  ans <- tibble()
  for (dd in list(...)) {
    
    if (! is.null(subsets)) {
      dd <- dd %>%
        filter(trait %in% subsets)
    }
    
    # total number of genes
    # total <- length(unique(dd$trait))
    
    dd <- dd %>%
      filter(!is.na(snp))
    
    # total number of genes that have CS
    cs_num <- length(unique(dd$trait))
    
    # average CS size
    tmp <- dd %>%
      group_by(trait, CSIndex) %>%
      summarize(n = n()) %>%
      ungroup() %>%
      summarize(n = mean(n))
    
    cs_size <- tmp$n
    
    # average pip
    tmp <- dd %>%
      group_by(trait, CSIndex) %>%
      summarize(n1 = mean(pip)) %>%
      ungroup() %>%
      summarize(n1 = mean(n1))
    
    avg_pip <- tmp$n1
    
    # average pip
    tmp <- dd %>%
      group_by(trait, CSIndex) %>%
      summarize(n1 = sum(any(pip>0.95))) %>%
      ungroup() %>%
      summarize(n1 = mean(n1))
    
    g09 <- tmp$n1
    
    ans <- bind_rows(ans, tibble(
      "total" = total,
      "cs_num" = cs_num,
      "cs_size" = cs_size,
      "avg_pip" = avg_pip,
      "g09" = g09))
    
  }
  return(ans)
}

basic_compare <- function(base, ..., subsets=NULL) {
  ans <- tibble()
  
  if (! is.null(subsets)) {
    base <- base %>%
      filter(trait %in% subsets)
  }
  
  new_base <- base %>%
    filter(!is.na(snp)) %>%
    group_by(trait, CSIndex) %>%
    summarize(cs_size = n(),
      avg_pip = mean(pip),
      g09 = sum(any(pip>0.95)))
  
  ct <- 0
  for (dd in list(...)) {
    ct <- ct + 1
    dd <- dd %>%
      filter(!is.na(snp))
    
    if (! is.null(subsets)) {
      dd <- dd %>%
        filter(trait %in% subsets)
    }
    
    new_dd <- dd %>%
      filter(!is.na(snp)) %>%
      group_by(trait, CSIndex) %>%
      summarize(cs_size = n(),
        avg_pip = mean(pip),
        g09 = sum(any(pip>0.95)))
    
    tmp_res <- new_base %>%
      inner_join(new_dd, by = c("trait", "CSIndex")) %>%
      group_by(trait) %>%
      summarize(cs_size.x = mean(cs_size.x),
        cs_size.y = mean(cs_size.y),
        avg_pip.x = mean(avg_pip.x),
        avg_pip.y = mean(avg_pip.y),
        g09.x = mean(g09.x),
        g09.y = mean(g09.y))
    
    c1 <- tidy(t.test(tmp_res$cs_size.x,
      tmp_res$cs_size.y)) %>%
      mutate(n = nrow(tmp_res),
        metric = "cs_size")
    
    c2 <- tidy(t.test(tmp_res$avg_pip.x,
      tmp_res$avg_pip.y)) %>%
      mutate(n = nrow(tmp_res),
        metric = "avg_pip")
    
    c3 <- tidy(t.test(tmp_res$g09.x,
      tmp_res$g09.y)) %>%
      mutate(n = nrow(tmp_res),
        metric = "g09")
    
    
    ans <- ans %>%
      bind_rows(bind_rows(c1, c2, c3) %>%
          mutate(type = ct) %>%
          select(estimate, p.value, n, metric, type))
  }
  
  return(ans)
}

prepare_corr <- function(df, sig, n_pop) {
  df_wk <- df[c("trait", "CSIndex", colnames(df)[grepl("corr", colnames(df))])] %>%
    pivot_longer(c(-trait, -CSIndex))
  
  single_sig <- sig %>%
    filter(p_value < 0.05)
  
  res <- tibble()
  
  res <- res %>%
    bind_rows(df_wk %>%
        group_by(name, CSIndex) %>%
        summarize(corr = mean(value),
          se = sd(value)/sqrt(n())) %>%
        filter(CSIndex == 1) %>%
        select(-CSIndex) %>%
        mutate(type = "1"))
  
  # simple sig
  res <- res %>%
    bind_rows(df_wk %>%
        filter(trait %in% unique(single_sig$trait)) %>%
        group_by(name, CSIndex) %>%
        summarize(corr = mean(value),
          se = sd(value)/sqrt(n())) %>%
        filter(CSIndex == 1) %>%
        select(-CSIndex) %>%
        mutate(type = "2"))
  
  # both sig
  new_corr <- tibble()
  for (idx in 1:(n_pop-1)) {
    for (jdx in (idx+1):n_pop) {
      df_tmp <- single_sig %>%
        filter(ancestry %in% c(idx, jdx)) %>%
        group_by(trait) %>%
        summarize(n = n()) %>%
        filter(n == 2)
      
      ans <- df_wk %>%
        filter(name %in% glue("ancestry{idx}_ancestry{jdx}_est_corr")) %>%
        filter(trait %in% unique(df_tmp$trait)) %>%
        group_by(name, CSIndex) %>%
        summarize(corr = mean(value),
          se = sd(value)/sqrt(n())) %>%
        filter(CSIndex == 1)
      
      new_corr <- new_corr %>%
        bind_rows(tibble(name = glue("ancestry{idx}_ancestry{jdx}_est_corr"),
          corr = ans$corr,
          se = ans$se))
    }
  }
  
  res <- res %>%
    bind_rows(new_corr %>%
        mutate(type = "3"))
  
  return(res)
}


