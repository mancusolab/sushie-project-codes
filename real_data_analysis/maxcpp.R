library(tidyverse)
library(glue)
library(broom)
library(stringr)
library(rmeta)

source("./utils.R")


# zenodo folder
ldsc_annotation_folder <- "~/Downloads/sushie_real_data_results/ldsc_annotation/"
ldsc_output_folder <- "~/Downloads/sushie_real_data_results/ldsc_output/"

traits <- c("EA_BAS", "EA_LYM", "EA_MON",
  "EA_WBC", "EA_NEU", "EA_EOS")

studies <- c("genoa", "mesa_rnaseq", "mesa_proteins")

types <- c("alpha")

methods <- c("sushie", "indep", "meta", "mega", "susiex", "mesusie")

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

all_res <- tibble()
for (study_name in studies) {
  anno_file <- read_tsv(glue("{ldsc_annotation_folder}/{study_name}.tsv.gz"))
  for (trait_name in traits) {
    for (method_name in methods) {
      for (type_name in types) {
        tmp_df <- read_tsv(glue("{ldsc_output_folder}/{trait_name}.{study_name}_{method_name}_{type_name}.bed.baselineLD.results")) %>%
          filter(Category %in% "L2_1") %>%
          select(Category, Prop._SNPs, Coefficient, Coefficient_std_error, `Coefficient_z-score`)
        colnames(tmp_df) <- c("cate", "p", "tau", "tau_se", "tau_z")
        tmp_log <- read_lines(glue("{ldsc_output_folder}/{trait_name}.{study_name}_{method_name}_{type_name}.bed.baselineLD.log")) 
        
        h2g <- as.numeric(gsub(" (.+)", "", gsub("Total Observed scale h2\\: ", "",
          tmp_log[grepl("scale h2", tmp_log)])))
        M <- as.numeric(str_extract(
          str_extract(tmp_log[grepl("SNPs with chi\\^2", tmp_log)],
            "\\(\\d+ SNPs remain\\)"), "\\d+"))
          anno_sd <- sd(
            c(anno_file[glue("{method_name}_{type_name}")][[1]],
              rep(0, M-nrow(anno_file))))
       
        
        all_res <- all_res %>%
          bind_rows(
            tmp_df %>%
              mutate(study = study_name,
                trait = trait_name,
                method = method_name,
                type = type_name,
                h2g = h2g,
                M = M,
                anno_sd = anno_sd)
          )
      }
    }
  }
}

summ_res <- tibble()
for (type_name in types) {
  for (method_name in methods) {
    tmp_all_res <- all_res %>%
      filter(type == type_name) %>%
      filter(method == method_name) %>%
      mutate(coef = anno_sd * M / h2g,
        tau_star = tau * coef,
        tau_star_se = tau_se * coef,
        tau_star_z = tau_star / tau_star_se)
    tmp_test <- meta.summaries(tmp_all_res$tau_star,
      tmp_all_res$tau_star_se, method="random")
    summ_res <- summ_res %>%
      bind_rows(
        tibble(type = type_name,
          method = method_name,
          tau_star = tmp_test$summary,
          tau_star_se = tmp_test$se.summary) %>%
          mutate(tau_star_z = tau_star / tau_star_se,
            uppci = tau_star + 1.96 * tau_star_se,
            lowci = tau_star - 1.96 * tau_star_se)
      )
  }
}

2*pnorm(6.73, lower.tail = F)

summ_res <- summ_res %>%
  mutate(method = factor(method,
    levels = c("sushie", "indep", "meta", "mega", "susiex", "mesusie"),
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE")))

ggplot(summ_res,
  aes(x = method, y = tau_star, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  geom_errorbar(aes(ymin = lowci, ymax = uppci),
    width = 0.1, position=position_dodge(width=0.5)) +
  geom_hline(yintercept = 0.448, linetype = "dashed") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 8, face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    axis.title.x=element_blank(),
    title = element_text(size = 10, face="bold"),
    axis.text=element_text(size = 8, face="bold"),
    text=element_text(size = 8))

# ggsave("./plots/s38.png", width = p_width-2, height = p_height+1)

enrich_comp <- summ_res %>%
  filter(method != "sushie") %>%
  ungroup() %>%
  mutate(sushie_tau_star = 0.448,
    sushie_tau_star_se = 0.0666,
    comp_z = (sushie_tau_star - tau_star) /
      sqrt(sushie_tau_star_se^2 + tau_star_se^2),
    p.value = 2*pnorm(abs(comp_z), lower.tail = F)) %>%
  select(method, tau_star, tau_star_se, sushie_tau_star,
    sushie_tau_star_se, comp_z, p.value)

