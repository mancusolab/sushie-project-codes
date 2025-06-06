library(tidyverse)
library(glue)
library(ggpubr)
library(broom)
library(RColorBrewer)

# to replicate our figures you need to download the data from the zenodo link
# and point it to simulation data paht
sim_data_path <- "~/Downloads/sushie_sim_data_results/"

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

method_colors <-c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a",
  "SuSiEx" = "#66a61e", "MESuSiE" = "#e6ab02", "XMAP" = "#a6761d", "XMAP-IND" = "#666666")

twas_colors <- c(method_colors,
  "LASSO" = "#fb8072", "Elastic Net" = "#b3de69", "gBLUP" = "#fccde5")

# 2 pop general data

load("./data/df_2pop.RData")

df_m1 <- df_2pop %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name, type) %>%
  filter(!is.na(value)) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == L2 & L3 == 0 & L1==2 &
      h2g %in% "0.05:0.05" & 
      rho %in% c("0.8")) %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

tmp_p1 <- filter(df_m1, type == "PIP")

small_p1 <- ggplot(tmp_p1,
  aes(x = factor(N), y = mvalue, color= factor(name))) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_color_manual(values = method_colors) +
  ylab("PIP of molQTLs") +
  xlab("molQTL Sample Size") +
  theme_sim()

tmp_p2 <- filter(df_m1, type == "CS")

small_p2 <- ggplot(tmp_p2,
  aes(x = factor(N), y = mvalue, color= factor(name))) +
  geom_point(size = point_size, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_color_manual(values =  method_colors) +
  scale_y_continuous(labels=scaleFUN) +
  xlab("molQTL Sample Size") +
  ylab("CS size") +
  theme_sim()

tmp_p3 <- filter(df_m1, type == "Calibration")

small_p3 <- ggplot(tmp_p3,
  aes(x = factor(N), y = mvalue, color= factor(name))) +
  geom_point(size = point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = method_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_y_continuous(labels=scaleFUN, breaks=c(seq(from=0,to=1,by=0.25), 0.9)) +
  ylab("Freq. of molQTLs in CS") +
  xlab("molQTL Sample Size") +
  theme_sim()

df_rho <- read_tsv(glue("{sim_data_path}/sushie_2pop_rho_all.tsv.gz"))  %>%
  filter(Lidx==1) %>%
  filter(method %in% c("sushie", "mesusie", "xmap", "xmap_ind")) %>%
  filter(N %in% "400:400" & L2 == 2 & L1 ==2 & L3==0 & h2g %in% "0.05:0.05") %>%
  group_by(method, sim, N, L1, L2, L3, h2g, rho) %>%
  summarize(mrho = mean(est_rho),
    se = (sd(est_rho) / sqrt(n())),
    n = n()) %>%
  ungroup() %>%
  select(sim, rho, mrho, name=method, se, n) %>%
  mutate(name = factor(name,
    levels = c("sushie", "mesusie", "xmap", "xmap_ind"),
    labels = c("SuShiE", "MESuSiE", "XMAP", "XMAP-IND")))

small_p4 <- ggplot(df_rho ,
  aes(x = factor(rho), y = mrho, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mrho -  1.96 *se, ymax = mrho +  1.96 *se),
    position=position_dodge(width=0.5), width = 0.2) +
  geom_hline(yintercept = c(0.01, 0.4, 0.8, 0.99),
    linetype="dashed",color="grey", alpha=0.5) +
  scale_y_continuous(labels = c("0.01", "0.4", "0.8", "0.99"),
    breaks = c(0.01, 0.4, 0.8, 0.99)) +
  scale_color_manual(values = method_colors) +
  ylab("Estimation") +
  xlab("True Effect Size Correlation") +
  theme_sim() 

load("./data/df_r2twas.RData")

df_r2 <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000
    & type == "R2") %>%
  group_by(sim, N, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(N = case_when(
      N == "200:200" ~ 200,
      N == "400:400" ~ 400,
      N == "600:600" ~ 600)) %>%
  filter(name %in% c("SuShiE", "SuShiE-Indep",
    "Meta-SuSiE", "SuSiE", "MESuSiE", 
    "XMAP", "XMAP-IND"))

small_p5 <- ggplot(df_r2,
  aes(x = factor(N), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous(labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("Training Sample Size") +
  theme_sim() 

df_twas <- df_r2twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & h2ge == 0.3/2000 & type == "TWAS") %>%
  group_by(sim, ngwas, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * sqrt(sd(value) / n())) %>%
  mutate(ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
      labels = c("100k", "200k", "300k"))) %>%
  filter(name %in% c("SuShiE", "SuShiE-Indep",
    "Meta-SuSiE", "SuSiE", "MESuSiE", 
    "XMAP", "XMAP-IND"))

small_p6 <- ggplot(df_twas,
  aes(x = factor(ngwas), y = mvalue, color = name)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("GWAS Sample Size") +
  theme_sim()

ggarrange(small_p1, small_p2, small_p3, small_p4, small_p5, small_p6,
  nrow = 3, ncol = 2, labels = c("a", "b", "c", "d", "e", "f"),
  common.legend = TRUE, legend="bottom", font.label = list(size=8))

# ggsave(filename = "./manuscript_plots/p2.pdf", width = p_width, height = 5.5)


