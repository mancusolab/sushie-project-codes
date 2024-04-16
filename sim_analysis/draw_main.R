library(tidyverse)
library(ggpubr)
library(broom)
library(RColorBrewer)

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

method_colors <- c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a")

twas_colors <- c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#a6cee3", "SuSiE" = "#e7298a",
  "LASSO" = "#66a61e", "Elastic Net" = "#e6ab02", "gBLUP" = "#a6761d")

pop2_pip <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_pip.tsv.gz")

ddPIP <- pop2_pip %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "PIP")

cp_pip <- pop2_pip %>%
  filter(h2g %in% "0.05:0.05" & L1 ==2) %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + h2g + rho, cp_pip)) %>%
  mutate(p.value = p.value/2)

pop2_cs <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_cs.tsv.gz")

ddCS <- pop2_cs %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "CS")

ss2 <- pop2_cs %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

pop2_prop <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_prop.tsv.gz")

ddPROP <- pop2_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  mutate(value = value/L2) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))

main_tt <- bind_rows(ddPIP, ddCS, ddPROP) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie"))) %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800))

tmp_p1 <- filter(main_tt, group == "PIP")

small_p1 <- ggplot(tmp_p1,
  aes(x = factor(N), y = mvalue, color= factor(name))) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
    position=position_dodge(width=0.5), width = 0.2) +
  scale_color_manual(values = method_colors) +
  ylab("PIP of molQTLs") +
  xlab("molQTL Sample Size") +
  theme_sim()

tmp_p2 <- filter(main_tt, group == "CS")

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

tmp_p3 <- filter(main_tt, group == "prec")

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


tmp_rho <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_rho.tsv.gz")

df_rho1 <- tmp_rho %>%
  filter(N %in% "400:400" & L2 == 2 & L1 ==2 & L3==0 & h2g %in% "0.05:0.05")

df_rho1 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

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

small_p4 <- ggplot(rho_all1,
  aes(x = rho, y = mrho, color = type)) +
  geom_abline(slope=1, intercept = 0) +
  geom_point() +
  geom_errorbar(aes(ymin = mrho - se, ymax = mrho + se),
    width = 0.02) +
  scale_y_continuous( labels=scaleFUN) +
  scale_x_continuous(breaks = c(0.01, 0.4, 0.8, 0.99)) +
  scale_color_brewer(palette = "Paired") +
  ylab("Est. Effect Size Correlation") +
  xlab("True Effect Size Correlation") +
  theme_sim() 


tmp_r2 <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_pred3_r2.tsv.gz")

df_r2 <- tmp_r2 %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  # filter(name == "ancestry2_weight2") %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, method) %>%
  summarize(adj_r2 = mean(value),
    adj_se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(method = factor(method,
    levels = c("sushie", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    N = case_when(
      N == "200:200" ~ 200,
      N == "400:400" ~ 400,
      N == "600:600" ~ 600,
      N == "800:800" ~ 800))

small_p5 <- ggplot(df_r2,
  aes(x = factor(N), y = adj_r2, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = adj_r2 - adj_se, ymax = adj_r2 + adj_se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("Training Sample Size") +
  theme_sim() 


cp_r2 <- tmp_r2 %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  filter(!is.na(value)) %>%
  mutate(method = factor(method,
    levels = c("sushie", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")))

tidy(lm(value ~ method + name + N, cp_r2)) %>%
  mutate(p.value = p.value/2)

tmp_twas <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_pred3_twas.tsv.gz")

df_twas <- tmp_twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & h2ge == 0.3/2000) %>%
  pivot_longer(cols=c(sushie, indep, meta, susie, enet, lasso, ridge)) %>%
  mutate(pval = 2 * pnorm(abs(value), lower.tail = FALSE)) %>%
  group_by(1:n()) %>%
  mutate(adjpval = min(1, pval*1000)) %>%
  # filter(ancestry == 2) %>%
  group_by(sim, h2ge, ngwas, name) %>%
  summarize(adj_power = mean(adjpval < 0.05),
    adj_se = 1.96 * sqrt(adj_power * (1 - adj_power) / n())) %>%
  filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name, levels = c("sushie", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
      labels = c("100k", "200k", "300k"))) %>%
  rename(method = name) %>%
  pivot_longer(cols = c(adj_power, adj_se)) %>%
  separate(name, sep = "_", into = c("type", "name")) %>%
  pivot_wider(values_from = value, names_from = "name")

small_p6 <- ggplot(df_twas,
  aes(x = factor(ngwas), y = power, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = power - se, ymax = power + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("GWAS Sample Size") +
  theme_sim()

cp_twas <- tmp_twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & h2ge == 0.3/2000) %>%
  pivot_longer(cols=c(sushie, indep, meta, susie, enet, lasso, ridge)) %>%
  mutate(value = value**2) %>%
  filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name, levels = c("sushie", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
      labels = c("100k", "200k", "300k")))

tidy(lm(value ~ name+ ngwas +factor(ancestry), cp_twas)) %>%
  mutate(p.value = p.value/2)

library(cowplot)
legend1 <- get_legend(small_p1)
legend2 <- get_legend(small_p4)
legend3 <- get_legend(small_p5)

prow1 <- plot_grid(small_p1 + theme(legend.position="none"),
  small_p2 + theme(legend.position="none"),
  small_p3 + theme(legend.position="none"),
  nrow = 1, align = "h",
  labels = c("A", "B", "C"), label_size=10)

prow2 <- plot_grid(small_p4 + theme(legend.position="none"),
  small_p5 + theme(legend.position="none"),
  small_p6 + theme(legend.position="none"),
  nrow = 1, align = "h",
  labels = c("D", "E", "F"), label_size= 10)

leg2 <- plot_grid(legend2, legend3, rel_widths = c(1, 2))

plot_grid(prow1, legend1, prow2, leg2, nrow=4, rel_heights = c(10,1,10,1))

ggsave(filename = "./manuscript_plots/main/p1.png",
  width = p_width, height = 5)




