library(tidyverse)
library(glue)
library(broom)
library(ggpubr)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Downloads/sushie_sim_data_results"

# 2 pop general data
pp_auprc <- read_tsv(glue("{sim_data_path}/auprc_data.tsv")) %>%
  mutate(name = factor(name,
    levels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

df_cali <- read_tsv(glue("{sim_data_path}/calibration_data.tsv")) %>%
  mutate(name = factor(name,
    levels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

ffont <- "sans"
fontsize <- 7
legend_fontsize <- 7
errbar_width <- 0.5
scaleFUN <- function(x) sprintf("%.2f", x)

p_width <- 7.08
p_height <- 2.5
point_size <- 1.5

theme_sim <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = legend_fontsize),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold", size=fontsize),
    axis.text=element_text(size = fontsize),
    legend.position = "bottom")
}

# auprc
df_tmp <- pp_auprc %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8" & h2g == "0.05:0.05") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

p1 <- ggplot(df_tmp, aes(x = Recall, y = Precision, color = name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~N, scales = "free") +
  scale_color_manual(values = method_colors) +
  theme_sim()

df_tmp <- pp_auprc %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2") &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8") %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2))

p2 <- ggplot(df_tmp, aes(x = Recall, y = Precision, color = name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~h2g, scales = "free") +
  scale_color_manual(values = method_colors) +
  theme_sim()

df_tmp <- pp_auprc %>%
  filter(h2g %in%  "0.05:0.05" &
      N %in%  "400:400" &
      L1 == L2  & L3 == 0 & rho == "0.8")

p3 <- ggplot(df_tmp, aes(x = Recall, y = Precision, color = name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~L2, scales = "free") +
  scale_color_manual(values = method_colors) +
  theme_sim()

df_tmp <- pp_auprc %>%
  filter(h2g %in% "0.05:0.05" &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 )

p4 <- ggplot(df_tmp, aes(x = Recall, y = Precision, color = name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~rho, scales = "free") +
  scale_color_manual(values = method_colors) +
  theme_sim()

ggarrange(p1, p2, p3 , p4, common.legend = TRUE, legend = "bottom",
  labels = c("eQTL Sample Size", "cis-SNP Heritability",
    "Number of cis-molQTLs", "Effect size correlation"),
  font.label = list(size = 8))

# ggsave("./manuscript_plots/additional/s1.png", width = 12, height = 8)

#  calibration
custom_labeller <- function(variable, value) {
  if (variable == "N") {
    return(as.character(value))  # Show only "N"
  } else {
    return("")  # Hide "name"
  }
}

df_tmp <- df_cali %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8" & h2g == "0.05:0.05")  %>%
  mutate(N = case_when(
    N == "200:200" ~ "N=200",
    N == "400:400" ~ "N=400",
    N == "600:600" ~ "N=600",
    N == "800:800" ~ "N=800"))

ggplot(df_tmp , aes(x = m_val, y = obs, color = name)) +
  geom_point() +
  geom_errorbar(aes(ymin = obs - 1.96*obs_se, ymax = obs + 1.96*obs_se), width = 0.05,) +
  geom_line() +
  facet_wrap(name~N, ncol = 4, labeller = labeller(N = label_value, name = ~"")) +
  scale_color_manual(values = method_colors) +
  geom_abline(intercept = 0, slope = 1, color="red",size=1, linetype="dashed") +
  xlab("Mean PIP") +
  ylab("Observed") +
  theme_sim()

# ggsave("./manuscript_plots/additional/s2.png", width = 10, height = 12)

df_tmp <- df_cali %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2") &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8") %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ "h2g=0.01",
    h2g == "0.05:0.05" ~ "h2g=0.05",
    h2g == "0.1:0.1" ~ "h2g=0.1",
    h2g == "0.2:0.2" ~ "h2g=0.2"))

ggplot(df_tmp , aes(x = m_val, y = obs, color = name)) +
  geom_point() +
  geom_errorbar(aes(ymin = obs - 1.96*obs_se, ymax = obs + 1.96*obs_se), width = 0.05,) +
  geom_line() +
  facet_wrap(name~h2g, ncol = 4, labeller = labeller(N = label_value, name = ~"")) +
  scale_color_manual(values = method_colors) +
  geom_abline(intercept = 0, slope = 1, color="red",size=1, linetype="dashed") +
  xlab("Mean PIP") +
  ylab("Observed") +
  theme_sim()

# ggsave("./manuscript_plots/additional/s3.png", width = 10, height = 12)

df_tmp <- df_cali %>%
  filter(h2g %in%  "0.05:0.05" &
      N %in%  "400:400" &
      L1 == L2  & L3 == 0 & rho == "0.8") %>%
  mutate(L2 = factor(case_when(
    L2 == 1 ~ "# of QTL: 1",
    L2 == 2 ~ "# of QTL: 2",
    L2 == 3 ~ "# of QTL 3"),
    levels = c("# of QTL: 1", "# of QTL: 2", "# of QTL 3")))

ggplot(df_tmp , aes(x = m_val, y = obs, color = name)) +
  geom_point() +
  geom_errorbar(aes(ymin = obs - 1.96*obs_se, ymax = obs + 1.96*obs_se), width = 0.05,) +
  geom_line() +
  facet_wrap(name~L2, ncol = 3, labeller = labeller(N = label_value, name = ~"")) +
  scale_color_manual(values = method_colors) +
  geom_abline(intercept = 0, slope = 1, color="red",size=1, linetype="dashed") +
  xlab("Mean PIP") +
  ylab("Observed") +
  theme_sim()

# ggsave("./manuscript_plots/additional/s4.png", width = 10, height = 12)

df_tmp <- df_cali %>%
  filter(h2g %in% "0.05:0.05" &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 ) %>%
  mutate(rho = factor(case_when(
    rho == 0 ~ "rho=0",
    rho == 0.01 ~ "rho=0.01",
    rho == 0.4 ~ "rho=0.4",
    rho == 0.8 ~ "rho=0.8",
    rho == 0.99 ~ "rho=0.99"),
    levels = c("rho=0", "rho=0.01", "rho=0.4", "rho=0.8", "rho=0.99")))

ggplot(df_tmp , aes(x = m_val, y = obs, color = name)) +
  geom_point() +
  geom_errorbar(aes(ymin = obs - 1.96*obs_se, ymax = obs + 1.96*obs_se), width = 0.05,) +
  geom_line() +
  facet_wrap(name~rho, ncol = 5, labeller = labeller(N = label_value, name = ~"")) +
  scale_color_manual(values = method_colors) +
  geom_abline(intercept = 0, slope = 1, color="red",size=1, linetype="dashed") +
  xlab("Mean PIP") +
  ylab("Observed") +
  theme_sim()

# ggsave("./manuscript_plots/additional/s5.png", width = 10, height = 12)

# cs vs generative L
df_cs_pop2_tmp <- read_tsv(glue("{sim_data_path}/sushie_2pop_cs_all.tsv.gz"))

df_num_cs <- df_cs_pop2_tmp %>%
  filter(L1 == L2) %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
  crossing(method = unique(df_cs_pop2_tmp$method)) %>%
  filter(method != "sushie_ss") %>%
  left_join(df_cs_pop2_tmp %>%
      filter(L1 == L2) %>%
      filter(method != "sushie_ss") %>%
      select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, SNPIndex_1based, causal,
        name = method) %>%
      filter(!is.na(SNPIndex_1based)) %>%
      distinct(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, name) %>%
      group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
      summarize(value = n()) %>%
      rename(method = name) %>%
      ungroup()) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, method) %>%
  summarize(mvalue = mean(value),
    se = sd(value)/sqrt(n())) %>%
  mutate(method= factor(method,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND"))) 

df_tmp <- df_num_cs %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8" & h2g == "0.05:0.05") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

p1 <- ggplot(df_tmp , aes(x = factor(N), y = mvalue, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - 1.96*se, ymax = mvalue + 1.96*se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Num. of CSs output") +
  ylim(0, 2.5) +
  xlab("cis-molQTL Sample Size with 2 causal QTLs") +
  theme_sim()

df_tmp <- df_num_cs %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2") &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8") %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2))

p2 <- ggplot(df_tmp , aes(x = factor(h2g), y = mvalue, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - 1.96*se, ymax = mvalue + 1.96*se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Num. of CSs output") +
  ylim(0, 2.5) +
  xlab("cis-SNP Heritability with 2 causal QTLs") +
  theme_sim()

df_tmp <- df_num_cs %>%
  filter(h2g %in% "0.05:0.05" &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 )

p3 <- ggplot(df_tmp , aes(x = factor(rho), y = mvalue, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - 1.96*se, ymax = mvalue + 1.96*se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Num. of CSs output") +
  ylim(0, 2.5) +
  xlab("Effect Size Correlation with 2 causal QTLs") +
  theme_sim()

df_tmp <- df_num_cs %>%
  filter(h2g %in%  "0.05:0.05" &
      N %in%  "400:400" &
      L1 == L2  & L3 == 0 & rho == "0.8")

p4 <- ggplot(df_tmp , aes(x = factor(L2), y = mvalue, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - 1.96*se, ymax = mvalue + 1.96*se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Num. of CSs output") +
  ylim(0, 2.5) +
  xlab("Number of cis-molQTL") +
  theme_sim()

ggarrange(p1, p2, p4,  p3,
  labels = c("A", "B", "C", "D"),
 common.legend = TRUE, legend = "bottom")

# ggsave("./manuscript_plots/additional/s6.png", width = 6, height = 4)

df_num_cs <- df_cs_pop2_tmp %>%
  filter(L1 == L2) %>%
  filter(L1 == 2) %>%
  filter(method != "sushie_ss") %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
  crossing(method = unique(df_cs_pop2_tmp$method)) %>%
  filter(method != "sushie_ss") %>%
  left_join(df_cs_pop2_tmp %>%
      filter(L1 == L2) %>%
      filter(method != "sushie_ss") %>%
      select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex,
        SNPIndex_1based, causal, name = method) %>%
      filter(!is.na(SNPIndex_1based)) %>%
      distinct(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, name) %>%
      group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
      summarize(value = n()) %>%
      rename(method = name) %>%
      ungroup()) %>%
  mutate(value = ifelse(is.na(value), 0, value))

df_res <- tibble()
for (method_name in c("mesusie", "indep", "meta", "susie", "susiex",
  "xmap", "xmap_ind")) {
  df_tmp <- df_num_cs %>%
    filter(method %in% c("sushie", method_name)) %>%
    mutate(method = factor(method,
      levels = c("sushie", method_name)))
  df_res <- df_res %>%
    bind_rows(
      tidy(lm(value ~ method + N + L1 +h2g +rho, df_tmp)) %>%
        filter(grepl("method", term))  
    )
}

df_res

df_res %>%
  filter(!grepl("susiex", term)) %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

# misspecified
df_tmp <- df_cs_pop2_tmp %>%
  filter(L2 >= L1 & L1 == 2) %>%
  filter(h2g %in% "0.05:0.05" & rho == 0.8 & N == "400:400" & L3 == 0) %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
  crossing(method = unique(df_cs_pop2_tmp$method)) %>%
  filter(method != "sushie_ss") %>%
  left_join(df_cs_pop2_tmp %>%
      filter(method != "sushie_ss") %>%
      select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex,
        SNPIndex_1based, causal, name = method) %>%
      filter(!is.na(SNPIndex_1based)) %>%
      distinct(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, name) %>%
      group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
      summarize(value = n()) %>%
      rename(method = name) %>%
      ungroup()) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, method) %>%
  summarize(mvalue = mean(value),
    se = sd(value)/sqrt(n())) %>%
  mutate(method= factor(method,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

ggplot(df_tmp , aes(x = factor(L2), y = mvalue, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - 1.96*se, ymax = mvalue + 1.96*se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Num. of CSs output") +
  ylim(0, 2.5) +
  xlab("Number of Single Effects Specified in Inference (L) with 2 causals") +
  theme_sim()

# ggsave("./manuscript_plots/additional/s7.png", width = 4, height = 3.5)

df_num_cs <- df_cs_pop2_tmp %>%
  filter(L1 <= L2) %>%
  filter(L1 == 2) %>%
  filter(h2g %in% "0.05:0.05" & rho == 0.8 & N == "400:400" & L3 == 0) %>%
  filter(method != "sushie_ss") %>%
  distinct(sim, locus, N, L1, L2, L3, h2g, rho) %>%
  crossing(method = unique(df_cs_pop2_tmp$method)) %>%
  filter(method != "sushie_ss") %>%
  left_join(df_cs_pop2_tmp %>%
      filter(L1 < L2) %>%
      filter(L1 == 2) %>%
      filter(method != "sushie_ss") %>%
      select(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, SNPIndex_1based, causal,
        name = method) %>%
      filter(!is.na(SNPIndex_1based)) %>%
      distinct(sim, locus, N, L1, L2, L3, h2g, rho, CSIndex, name) %>%
      group_by(sim, locus, N, L1, L2, L3, h2g, rho, name) %>%
      summarize(value = n()) %>%
      rename(method = name) %>%
      ungroup()) %>%
  mutate(value = ifelse(is.na(value), 0, value))

df_res <- tibble()
for (method_name in c("mesusie", "indep", "meta", "susie", "susiex",
  "xmap", "xmap_ind")) {
  df_tmp <- df_num_cs %>%
    # filter(value <=2) %>%
    filter(method %in% c("sushie", method_name)) %>%
    mutate(method = factor(method,
      levels = c("sushie", method_name)))
  df_res <- df_res %>%
    bind_rows(
      tidy(lm(value ~ method + L2, df_tmp)) %>%
        filter(grepl("method", term))  
    )
}

df_res

df_res %>%
  filter(!grepl("susiex", term)) %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# FDR-based power
tmp_df_pip_fdr <- read_tsv(glue("{sim_data_path}/fdr_pip_data.tsv"))

df_pip_fdr <- tmp_df_pip_fdr %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(Power_mean = mean(value),
    Power_se = sd(value)/sqrt(n())) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
      "XMAP", "XMAP-IND")))

df_tmp <- df_pip_fdr %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8" & h2g == "0.05:0.05") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

p1 <- ggplot(df_tmp, aes(x = factor(N), y = Power_mean, color = name)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Power") +
  ylim(0, 1) +
  xlab("cis-molQTL Sample Size") +
  theme_sim()

df_tmp <- df_pip_fdr %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2") &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8") %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2))

p2 <- ggplot(df_tmp, aes(x = factor(h2g), y = Power_mean, color = name)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Power") +
  ylim(0, 1) +
  xlab("cis-SNP Heritability") +
  theme_sim()

df_tmp <- df_pip_fdr %>%
  filter(h2g %in%  "0.05:0.05" &
      N %in%  "400:400" &
      L1 == L2  & L3 == 0 & rho == "0.8")

p3 <- ggplot(df_tmp, aes(x = factor(L2), y = Power_mean, color = name)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Power") +
  ylim(0, 1) +
  xlab("Number of cis-molQTL") +
  theme_sim()

df_tmp <- df_pip_fdr %>%
  filter(h2g %in% "0.05:0.05" &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 )

p4 <- ggplot(df_tmp, aes(x = factor(rho), y = Power_mean, color = name)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Power") +
  ylim(0, 1) +
  xlab("Effect Size Correlation") +
  theme_sim()

ggarrange(p1, p2, p3 , p4, common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D"), font.label = list(size = 8))

# ggsave("./manuscript_plots/additional/s8.png", width = 6, height = 4)

df_res <- tibble()
for (method_name in c("mesusie", "indep", "meta", "susie", "susiex",
  "xmap", "xmap_ind")) {
  df_tmp <- tmp_df_pip_fdr %>%
    filter(name %in% c("sushie", method_name)) %>%
    mutate(name = factor(name,
      levels = c("sushie", method_name)))
  
  df_res <- df_res %>%
    bind_rows(
      tidy(lm(value ~ name + N+L1+L2+L3+rho+h2g, df_tmp)) %>%
        filter(grepl("name", term))  
    )
}

df_res

df_res %>%
  filter(!grepl("susiex", term)) %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

tmp_df_cs_fdr <- read_tsv(glue("{sim_data_path}/fdr_cs_data.tsv"))

df_cs_fdr <- tmp_df_cs_fdr %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, method) %>%
  summarize(n =n(),
    Power_mean = mean(Power),
    Power_se = sd(Power)/sqrt(n())) %>%
  mutate(method= factor(method,
    levels = c("sushie", "indep", "meta", "susie",
      "susiex", "mesusie",  "xmap",  "xmap_ind"),
    labels = c("SuShiE", "SuShiE-Indep",
      "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE",
      "XMAP", "XMAP-IND")))
  
df_tmp <- df_cs_fdr %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8" & h2g == "0.05:0.05") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

p1 <- ggplot(df_tmp, aes(x = factor(N), y = Power_mean, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Coverage") +
  ylim(0,1) +
  xlab("cis-molQTL Sample Size") +
  theme_sim()

df_tmp <- df_cs_fdr %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2") &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 & rho == "0.8") %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2))

p2 <- ggplot(df_tmp, aes(x = factor(h2g), y = Power_mean, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Coverage") +
  ylim(0,1) +
  xlab("cis-SNP Heritability") +
  theme_sim()

df_tmp <- df_cs_fdr %>%
  filter(h2g %in%  "0.05:0.05" &
      N %in%  "400:400" &
      L1 == L2  & L3 == 0 & rho == "0.8")

p3 <- ggplot(df_tmp, aes(x = factor(L2), y = Power_mean, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Coverage") +
  ylim(0,1) +
  xlab("Number of cis-molQTL") +
  theme_sim()

df_tmp <- df_cs_fdr %>%
  filter(h2g %in% "0.05:0.05" &
      N %in%  "400:400" &
      L1 == 2 & L2 == 2 & L3 == 0 )

p4 <- ggplot(df_tmp, aes(x = factor(rho), y = Power_mean, color = method)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Power_mean - 1.96*Power_se, ymax = Power_mean + 1.96*Power_se),
    width = 0.05, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  ylab("Coverage") +
  ylim(0,1) +
  xlab("Effect Size Correlation") +
  theme_sim()

ggarrange(p1, p2, p3 , p4, common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D"), font.label = list(size = 8))

# ggsave("./manuscript_plots/additional/s9.png", width = 6, height = 4)

df_res <- tibble()
for (method_name in c("mesusie", "indep", "meta", "susie", "susiex",
  "xmap", "xmap_ind")) {
  df_tmp <- tmp_df_cs_fdr %>%
    filter(method %in% c("sushie", method_name)) %>%
    mutate(method = factor(method,
      levels = c("sushie", method_name)))
  
  df_res <- df_res %>%
    bind_rows(
      tidy(lm(Power ~ method + N+L1+L2+L3+rho+h2g, df_tmp)) %>%
        filter(grepl("method", term))  
    )
}

df_res

df_res %>%
  filter(!grepl("susiex", term)) %>%
  filter(!grepl("xmap", term)) %>%
  mutate(weight = 1/(std.error^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

