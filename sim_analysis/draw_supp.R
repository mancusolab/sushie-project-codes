library(tidyverse)
library(ggpubr)
library(broom)
library(RColorBrewer)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Documents/github/data/sushie_results/sim3"

ffont <- "sans"
fontsize <- 9
errbar_width <- 0.5
scaleFUN <- function(x) sprintf("%.2f", x)

p_width <- 7
p_height <- 2.5
point_size <- 1.5

theme_sim <- function() {
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.text = element_text(face="bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    title=element_text(face="bold"),
    text=element_text(size = fontsize))
}

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

twas_colors <- c(method_colors,
  "LASSO" = "#fb8072", "Elastic Net" = "#b3de69", "gBLUP" = "#fccde5")

pop3_colors <- c("One-Ancestry" = "#7fc97f", "Two-Ancestry" = "#beaed4",
  "Three-Ancestry" = "#fdc086")

pp <- function(tt1, chvar, xlab, colors_type = method_colors){
  envar = sym(chvar)
  
  p1tmp <- filter(tt1, type == "PIP")
  
  p1 <- ggplot(p1tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size=point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    ylab("PIP of molQTLs") +
    xlab(xlab) +
    theme_sim()
  
  p2tmp <- filter(tt1, type == "CS")
  
  p2 <- ggplot(p2tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    scale_y_continuous(labels=scaleFUN) +
    ylab("CS size") +
    xlab(xlab) +
    theme_sim()
  
  p3tmp <- filter(tt1, type == "Calibration")
  
  p3 <- ggplot(p3tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_y_continuous(labels=scaleFUN, breaks=c(seq(from=0,to=1,by=0.25), 0.9)) +
    ylab("Freq. of molQTLs in CS") +
    xlab(xlab) +
    theme_sim()
  
  p <- ggarrange(p1, p2, p3, labels = c("A", "B","C"), nrow = 1,
    font.label = list(size = 8), common.legend = TRUE, legend = "bottom")
  return(p)
}

# 3 population
load("./data/df_3pop.RData")

df_tmp <- df_3pop %>%
  group_by(name, type, N) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) 

s_p1 <- pp(df_tmp, "N", "molQTL Sample Size", colors_type = pop3_colors)

# ggsave("./manuscript_plots/supp/s1.png", width = p_width, height = p_height+0.5)
load("./data/df_2pop.RData")

df_2pop %>%
  filter(h2g == "0:0" & type == "CS") %>%
  distinct(name, locus) %>%
  group_by(name) %>%
  summarize(n = n()) %>%
  mutate(prop = n / 500,
    z = (prop - 0.01) / sqrt(0.01 * (1 - 0.01) / 500))

general_sim1 <- function(tt1, ylab, colors_type = method_colors) {
  p2tmp <- tt1 %>%
    filter(N == "400:400" & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
    filter(L1 == L2)
  
  p2 <- ggplot(p2tmp,
    aes(x = factor(L1), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    scale_y_continuous(labels=scaleFUN, limits = c(0, 1)) +
    ylab(ylab) +
    xlab("The number of molQTLs") +
    theme_sim()
  
  p3tmp <- tt1 %>%
    filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
    filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
    mutate(h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2))
  
  p3 <- ggplot(p3tmp,
    aes(x = factor(h2g), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_y_continuous(labels=scaleFUN, limits = c(0,1)) +
    ylab(ylab) +
    xlab("cis-SNP Heritability") +
    theme_sim()
  
  p4tmp <- tt1 %>%
    filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05")
  
  p4 <- ggplot(p4tmp,
    aes(x = factor(rho), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_y_continuous(labels=scaleFUN, limits = c(0,1)) +
    ylab(ylab) +
    xlab("Effect size corr. across ancestries") +
    theme_sim()
  
  p <- ggarrange(p2, p3, p4, labels = c("A", "B", "C"), nrow=1,
    font.label = list(size = 8), common.legend = TRUE, legend = "bottom")
  return(p)
}


# sample size, l1, h2g, rho
df_tmp <- df_2pop %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  # filter(rho == 0.8) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  filter(!is.na(value)) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

general_sim1(filter(df_tmp, type == "PIP"), "PIP of molQTLs")

# ggsave("./manuscript_plots/supp/s2.png", width = p_width+0.4, height = p_height+0.5)

general_sim2 <- function(tt1, ylab, colors_type = method_colors){
  p2tmp <- tt1 %>%
    filter(N == "400:400" & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
    filter(L1 == L2)
  
  p2 <- ggplot(p2tmp,
    aes(x = factor(L1), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    # scale_y_continuous(labels=scaleFUN, limits = c(0, 7)) +
    ylab(ylab) +
    xlab("The number of molQTLs") +
    theme_sim()
  
  p3tmp <- tt1 %>%
    filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
    filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
    mutate(h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2))
  
  p3 <- ggplot(p3tmp,
    aes(x = factor(h2g), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    # scale_y_continuous(labels=scaleFUN, limits = c(0, 7)) +
    ylab(ylab) +
    xlab("cis-SNP Heritability") +
    theme_sim()
  
  p4tmp <- tt1 %>%
    filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05")
  
  p4 <- ggplot(p4tmp,
    aes(x = factor(rho), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    # scale_y_continuous(labels=scaleFUN, limits = c(0, 7)) +
    ylab(ylab) +
    xlab("Effect size corre. across ancestries") +
    theme_sim()
  
  p <- ggarrange(p2, p3, p4, labels = c("A", "B", "C"), nrow=1,
    font.label = list(size = 8), common.legend = TRUE, legend = "bottom")
  return(p)
}

general_sim2(filter(df_tmp, type == "CS"), "Credible Set Size")

# ggsave("./manuscript_plots/supp/s3.png", width = p_width+0.4, height = p_height+0.5)

general_sim1(filter(df_tmp, type == "Calibration"), "Freq. of molQTLs in CS")

# ggsave("./manuscript_plots/supp/s4.png", width = p_width+0.4, height = p_height+0.5)


# 2 population

tt1 <- df_2pop  %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  mutate(N = case_when(
    N == "400:200" ~ 200,
    N == "400:400" ~ 400,
    N == "400:600" ~ 600,
    N == "400:800" ~ 800)) %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

pp(tt1, "N", "Sample Size for 2nd Ancestry")

# ggsave("./manuscript_plots/supp/s5.png", width = p_width, height = p_height+0.5)

# different h2g

tt2 <- df_2pop %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  mutate(h2g = case_when(
    h2g == "0.05:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.05:0.1" ~ 0.1,
    h2g == "0.05:0.2" ~ 0.2))  %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

pp(tt2, "h2g", "cis-SNP h2g for 2nd Ancestry")

# ggsave("./manuscript_plots/supp/s6.png", width = p_width, height = p_height+0.5)

# rho
tmp_rho <- read_tsv(glue("{sim_data_path}sushie_2pop_rho.tsv.gz"))

sushie_cs_pop2 <- read_tsv(glue("{sim_data_path}/sushie_2pop_cs.tsv.gz"))

all_sim1 <- sushie_cs_pop2 %>%
  filter(N %in% "400:400" & L1 == L2 & L1 ==2 & L3==0 & h2g %in% "0.05:0.05") %>%
  filter(!is.na(sushie)) %>%
  distinct(sim, locus)

df_rho1 <- tmp_rho %>%
  filter(method %in% "sushie") %>%
  inner_join(all_sim1, by = c("sim", "locus"))

rho_11 <- df_rho1 %>%
  group_by(sim, rho) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n()))) %>%
  mutate(type ="All CSs")

rho_12 <- df_rho1 %>%
  filter(Lidx==1) %>%
  group_by(sim, rho) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n()))) %>%
  mutate(type ="First CSs")

rho_13 <- df_rho1 %>%
  filter(Lidx==2) %>%
  group_by(sim, rho) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n()))) %>%
  mutate(type ="Second CSs")

rho_all1 <- bind_rows(rho_11, rho_12, rho_13)

p_rho1 <- ggplot(rho_all1,
  aes(x = factor(rho), y = mrho, color = type)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mrho - se, ymax = mrho + se),
    width = 0.02, position=position_dodge(width=0.5)) +
  scale_y_continuous(breaks = c(0.01, 0.4, 0.8, 0.99), labels=scaleFUN) +
  geom_hline(yintercept = c(0.01, 0.4, 0.8, 0.99),
    linetype="dashed",color="grey", alpha=0.5) +
  scale_color_brewer(palette = "Paired") +
  ylab("Est. Effect Size Correlation") +
  xlab("True Effect Size Correlation") +
  theme_sim() 

all_sim2 <- sushie_cs_pop2 %>%
  filter(N %in% c("400:400", "1200:1200", "2400:2400") & L1 == L2 & L1 ==2 & L3==0 & h2g %in% "0.05:0.05" & rho == 0.8) %>%
  filter(!is.na(sushie)) %>%
  distinct(sim, locus)

df_rho2 <- tmp_rho %>%
  filter(method %in% "sushie") %>%
  inner_join(all_sim2, by = c("sim", "locus"))

rho_21 <- df_rho2 %>%
  group_by(sim, N) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n()))) %>%
  mutate(type ="All CSs")

rho_22 <- df_rho2 %>%
  filter(Lidx==1) %>%
  group_by(sim, N, rho) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n()))) %>%
  mutate(type ="First CSs")

rho_23 <- df_rho2 %>%
  filter(Lidx==2) %>%
  group_by(sim, N, rho) %>%
  summarize(mrho = mean(est_rho),
    se = 1.96 * (sd(est_rho) / sqrt(n())),
    n = n()) %>%
  mutate(type ="Second CSs")

rho_all2 <- bind_rows(rho_21, rho_22, rho_23) %>%
  mutate(N = case_when(
    N == "400:400" ~ 400,
    N == "1200:1200" ~ 1200,
    N == "2400:2400" ~ 2400))

p_rho2 <- ggplot(rho_all2,
  aes(x = factor(N), y = mrho, color = type)) +
  geom_hline(yintercept = 0.8) +
  geom_point( position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mrho - se, ymax = mrho + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous(breaks = c(seq(0, 1, by=0.25), 0.8), labels=scaleFUN) +
  xlab("molQTL Sample Size") +
  ylab("Estimated Correlation") +
  scale_color_brewer(palette = "Paired") +
  theme(panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(face = "bold"),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    axis.title=element_text(face="bold"),
    text=element_text(size = 10))

ggarrange(p_rho1, p_rho2, ncol = 2, labels = c("A", "B"),
  common.legend = TRUE, legend = "bottom")
# ggsave("./manuscript_plots/supp/s7.png", width = 5, height = 3)

tt3 <- df_2pop %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05" & L3 != 3) %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

tt3 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

pp(tt3, "L3", "Number of AS molQTLs")
# ggsave("./manuscript_plots/supp/s8.png", width = 10, height = 4)

# fdr
load("./data/df_noshared.RData")

df_tmp <- df_noshared %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) 

tmp_noshared1 <- df_tmp %>%
  filter(L1 == 0 & L2 == 4 & L3 == 2 & h2g == "0.05:0.05" & rho == 0) %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

p_noshared1 <- ggplot(tmp_noshared1,
  aes(x = factor(N), y = mvalue, color= factor(name))) +
  geom_point(size = point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab("Freq. CS capture non-causal molQTLs") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), limits=c(0,1)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  xlab("molQTL Sample size") +
  theme_sim()

tmp_noshared2 <- df_tmp %>%
  filter(L1 == 0 & L2 == L3 *2  & h2g == "0.05:0.05" & rho == 0) %>%
  filter(N %in%  "400:400")

p_noshared2 <- ggplot(tmp_noshared2,
  aes(x = factor(L3), y = mvalue, color= factor(name))) +
  geom_point(size = point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab("Freq. CS capture non-causal molQTLs") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), limits=c(0,1))+
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  xlab("Number of AS molQTLs") +
  theme_sim()

tmp_noshared3 <- df_tmp %>%
  filter(L1 == 0 & L2 == 4 & L3 ==2  & N %in%  "400:400" & rho == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2))

p_noshared3 <- ggplot(tmp_noshared3,
  aes(x = factor(h2g), y = mvalue, color= factor(name))) +
  geom_point(size = point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  ylab("Freq. CS capture non-causal molQTLs") +
  scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1), limits=c(0,1))+
  xlab("cis-SNP Heritability") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_sim()

ggarrange(p_noshared1, p_noshared2, p_noshared3,
  labels = c("A", "B", "C"), nrow=1,
  font.label = list(size = 8), common.legend = TRUE, legend = "bottom")

# ggsave("./manuscript_plots/supp/s9.png", width = p_width, height = p_height+0.75)

pp_l2 <- function(tt1, chvar, xlab, colors_type = method_colors){
  envar = sym(chvar)
  
  p1tmp <- filter(tt1, type == "PIP")
  
  p1 <- ggplot(p1tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size=point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    scale_y_continuous(labels=scaleFUN, limits=c(0,1)) +
    ylab("PIP of molQTLs") +
    xlab(xlab) +
    theme_sim()
  
  p2tmp <- filter(tt1, type == "CS")
  
  p2 <- ggplot(p2tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    scale_y_continuous(labels=scaleFUN) +
    ylab("CS size") +
    xlab(xlab) +
    theme_sim()
  
  p3tmp <- filter(tt1, type == "Calibration")
  
  p3 <- ggplot(p3tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size = point_size, position=position_dodge(width=0.5)) +
    scale_color_manual(values = colors_type) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_y_continuous(labels=scaleFUN, limits=c(0,1)) +
    # geom_hline(yintercept = 0.9, linetype = "dashed") +
    ylab("Freq. of molQTLs in CS") +
    xlab(xlab) +
    theme_sim()
  
  p <- ggarrange(p1, p2, p3, labels = c("A", "B","C"), nrow = 1,
    font.label = list(size = 8), common.legend = TRUE, legend = "bottom")
  return(p)
}

tt4 <- df_2pop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

tt4 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

pp_l2(tt4, "L2", "Number of Inferred molQTLs")

# ggsave("./manuscript_plots/supp/s10.png", width = p_width, height = p_height+0.5)

# twas r2 
load("./data/df_r2twas.RData")
  
df_r2_n <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000 & type == "R2") %>%
  group_by(sim, N, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(N = case_when(
      N == "200:200" ~ 200,
      N == "400:400" ~ 400,
      N == "600:600" ~ 600,
      N == "800:800" ~ 800))

r2_n <- ggplot(df_r2_n,
  aes(x = factor(N), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("Training Sample Size") +
  theme_sim() 

df_r2_h2g <- df_r2twas %>%
  filter(N == "400:400" &
      ngwas == 2e5 & h2ge == 0.3/2000 & type == "R2") %>%
  group_by(sim, h2g, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2))

r2_h2g <- ggplot(df_r2_h2g,
  aes(x = factor(h2g), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("cis-SNP Heritability") +
  theme_sim() 

df_twas_n <- df_r2twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & h2ge == 0.3/2000 & type == "TWAS") %>%
  group_by(sim, ngwas, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sqrt(sd(value) / n())) %>%
  mutate(ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
    labels = c("100k", "200k", "300k")))

twas_n <- ggplot(df_twas_n,
  aes(x = factor(ngwas), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("GWAS Sample Size") +
  theme_sim()

df_twas_h2ge <- df_r2twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & ngwas == 200000 & type == "TWAS") %>%
  group_by(sim, h2ge, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sqrt(sd(value) / n())) %>%
  mutate(h2ge = factor(case_when(
      h2ge == unique(tmp_twas$h2ge)[2] ~ "6e-5",
      h2ge == 0.00015 ~ "1.5e-4",
      h2ge == 0.00030 ~ "3e-4",
      h2ge == 0.0006 ~ "6e-4",
    ), levels = c("6e-5",  "1.5e-4", "3e-4", "6e-4")))

twas_h2ge <- ggplot(df_twas_h2ge,
  aes(x = factor(h2ge), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("Prop. of trait heritability mediated by expression") +
  theme_sim()

ggarrange(r2_n, r2_h2g, twas_n, twas_h2ge,
  common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D"), font.label = list(size = 8))

# ggsave("./manuscript_plots/supp/s11.png", width = p_width, height = p_height+2)

