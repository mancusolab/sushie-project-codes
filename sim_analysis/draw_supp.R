library(tidyverse)
library(ggpubr)
library(broom)
library(RColorBrewer)

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
    text=element_text(size = fontsize))
}

method_colors <- c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#7570b3", "SuSiE" = "#e7298a")

twas_colors <- c("SuShiE" = "#1b9e77", "SuShiE-Indep" = "#d95f02",
  "Meta-SuSiE" = "#a6cee3", "SuSiE" = "#e7298a",
  "LASSO" = "#66a61e", "Elastic Net" = "#e6ab02", "gBLUP" = "#a6761d")

pop3_colors <- c("One-Ancestry" = "#7fc97f", "Two-Ancestry" = "#beaed4",
  "Three-Ancestry" = "#fdc086")


pp <- function(tt1, chvar, xlab, colors_type = method_colors){
  envar = sym(chvar)
  
  p1tmp <- filter(tt1, group == "PIP")
  
  p1 <- ggplot(p1tmp,
    aes(x = factor(!!envar), y = mvalue, color= factor(name))) +
    geom_point(size=point_size, position=position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin = mvalue - se, ymax = mvalue + se),
      position=position_dodge(width=0.5), width = 0.2) +
    scale_color_manual(values = colors_type) +
    ylab("PIP of molQTLs") +
    xlab(xlab) +
    theme_sim()
  
  p2tmp <- filter(tt1, group == "CS")
  
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
  
  p3tmp <- filter(tt1, group == "prec")
  
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
pop3_pip <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_3pop_pip.tsv.gz")

ddPIP3 <- pop3_pip %>%
  select(sushie1, sushie2, sushie3, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie1:sushie3) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "PIP")

pop3_cs <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_3pop_cs.tsv.gz")

ddCS3 <- pop3_cs %>%
  select(sushie1, sushie2, sushie3, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie1 != 0 & sushie2 != 0 & sushie3 != 0) %>%
  pivot_longer(cols = sushie1:sushie3) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "CS")

pop3_prop <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_3pop_prop.tsv.gz")

ddPROP3 <- pop3_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L1) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n())))


ttdd <- bind_rows(ddPIP3, ddCS3, ddPROP3) %>%
  mutate(name = factor(name,
    labels = c("One-Ancestry", "Two-Ancestry", "Three-Ancestry"),
    levels = c("sushie1", "sushie2", "sushie3")))

ttdd_1 <- ttdd %>%
  filter(h2g == 0.05 & L2 == L1 & L2 == 2)

ttdd_1 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

main_p5 <- pp(ttdd_1, "N", "molQTL Sample Size", colors_type = pop3_colors)

# ggsave("./manuscript_plots/supp/s1.png", width = p_width, height = p_height+0.5)

cp1 <- pop3_pip %>%
  select(sushie1, sushie2, sushie3, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie1:sushie3) %>%
  mutate(name = factor(name, levels = c("sushie3", "sushie1", "sushie2")))

tidy(lm(value ~ name + N , cp1)) %>%
  mutate(p.value = p.value/2)

cp2 <- pop3_cs %>%
  select(sushie1, sushie2, sushie3, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie1 != 0 & sushie2 != 0 & sushie3 != 0) %>%
  pivot_longer(cols = sushie1:sushie3) %>%
  mutate(name = factor(name, levels = c("sushie3", "sushie1", "sushie2")))

tidy(lm(value ~ name + N, cp2)) %>%
  mutate(p.value = p.value/2)

cp3 <- pop3_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L1) %>%
  mutate(name = factor(name, levels = c("sushie3", "sushie1", "sushie2")))

tidy(lm(value ~ name + N, cp3)) %>%
  mutate(p.value = p.value/2)

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
    scale_y_continuous(labels=scaleFUN, limits = c(0,1)) +
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

pop2_pip <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_pip.tsv.gz")

# sample size, l1, h2g, rho
df_pip <- pop2_pip %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie"))) 

general_sim1(df_pip, "PIP of molQTLs")

# ggsave("./manuscript_plots/supp/s2.png", width = p_width+0.4, height = p_height+0.5)

cp_pip <- pop2_pip %>%
  filter(N %in% "400:400") %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + h2g + rho, cp_pip)) %>%
  mutate(p.value = p.value/2)

all_pip <- pop2_pip %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + L2 + L3 + N + h2g + rho, all_pip)) %>%
  mutate(p.value = p.value/2) %>%
  filter(term %in% c("nameindep", "namemeta", "namesusie")) %>%
  select(estimate) %>%
  unlist() %>%
  as.numeric() %>%
  mean()

cp_pip2 <- pop2_pip %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
  filter(h2g %in% "0.05:0.05") %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + N, cp_pip2)) %>%
  mutate(p.value = p.value/2)

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

pop2_cs <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_cs.tsv.gz")

df_cs <- pop2_cs %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

general_sim2(df_cs, "Credible Set Size")

# ggsave("./manuscript_plots/supp/s3.png", width = p_width+0.4, height = p_height+0.5)

cp_cs <- pop2_cs %>%
  filter(N %in% "400:400") %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + h2g + rho, cp_cs)) %>%
  mutate(p.value = p.value/2)

all_cs <- pop2_cs %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name, levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + L2 + L3 + N + h2g + rho, all_cs)) %>%
  mutate(p.value = p.value/2) %>%
  filter(term %in% c("nameindep", "namemeta", "namesusie")) %>%
  select(estimate) %>%
  unlist() %>%
  as.numeric() %>%
  mean()

cp_cs2 <- pop2_cs %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
  filter(h2g %in% "0.05:0.05") %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + N, cp_cs2)) %>%
  mutate(p.value = p.value/2)

pop2_prop <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_prop.tsv.gz")

df_prop <- pop2_prop %>%
  filter(N %in% "400:400") %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  mutate(value = value/L2) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

general_sim1(df_prop, "Freq. of molQTLs in CS")

# ggsave("./manuscript_plots/supp/s4.png", width = p_width+0.4, height = p_height+0.5)

cp_prop <- pop2_prop %>%
  filter(N %in% "400:400") %>%
  filter(L1 == L2 & L3 == 0) %>%
  filter(h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + h2g + rho, cp_prop)) %>%
  mutate(p.value = p.value/2)

all_prop <- pop2_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L1 + L2 + L3 + N + h2g + rho, all_prop)) %>%
  mutate(p.value = p.value/2) %>%
  filter(term %in% c("nameSuShiE-Indep", "nameMeta-SuSiE", "nameMega-SuSiE")) %>%
  select(estimate) %>%
  unlist() %>%
  as.numeric() %>%
  mean()

cp_prop2 <- pop2_prop %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800")) %>%
  filter(L1 == L2 & L3 == 0 & L1 == 2) %>%
  filter(h2g %in% "0.05:0.05") %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))


tidy(lm(value ~ name + N, cp_prop2)) %>%
  mutate(p.value = p.value/2)


# 2 population
pop2_pip <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_pip.tsv.gz")

ddPIP2 <- pop2_pip %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "PIP")

pop2_cs <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_cs.tsv.gz")

ddCS2 <- pop2_cs %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  group_by(sim, N, L1, L2, L3, h2g, rho, name) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "CS")

pop2_prop <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_prop.tsv.gz")

ddPROP2 <- pop2_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L2) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "prec")


ttdd <- bind_rows(ddPIP2, ddCS2, ddPROP2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie"))) 


# sample size
tt1 <- ttdd %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  mutate(N = case_when(
    N == "400:200" ~ 200,
    N == "400:400" ~ 400,
    N == "400:600" ~ 600,
    N == "400:800" ~ 800))

pp(tt1, "N", "Sample Size for 2nd Ancestry")

# ggsave("./manuscript_plots/supp/s5.png", width = p_width, height = p_height+0.5)

cp_ss1 <- pop2_pip %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + N, cp_ss1)) %>%
  mutate(p.value = p.value/2)

cp_ss2 <- pop2_cs %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + N, cp_ss2)) %>%
  mutate(p.value = p.value/2)

cp_ss3 <- pop2_prop %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "Mega-SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + N, cp_ss3)) %>%
  mutate(p.value = p.value/2)

# different h2g
tt2 <- ttdd %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  mutate(h2g = case_when(
    h2g == "0.05:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.05:0.1" ~ 0.1,
    h2g == "0.05:0.2" ~ 0.2))

pp(tt2, "h2g", "cis-SNP h2g for 2nd Ancestry")

# ggsave("./manuscript_plots/supp/s6.png", width = p_width, height = p_height+0.5)

cp_ss1 <- pop2_pip %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + h2g, cp_ss1)) %>%
  mutate(p.value = p.value/2)

cp_ss2 <- pop2_cs %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + h2g, cp_ss2)) %>%
  mutate(p.value = p.value/2)

cp_ss3 <- pop2_prop %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + h2g, cp_ss3)) %>%
  mutate(p.value = p.value/2)

# rho
tmp_rho <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_2pop_rho.tsv.gz")


df_rho2 <- tmp_rho %>%
  filter(N %in% c("400:400", "1200:1200", "2400:2400") & L1 == L2 & L1 ==2 & L3==0 & h2g %in% "0.05:0.05" & rho == 0.8)

df_rho2 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n=n())

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

ggplot(rho_all2,
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
    text=element_text(size = 6))

# ggsave("./manuscript_plots/supp/s7.png", width = 3, height = 3)

# ancestry-specific
tt3 <- ttdd %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05" & L3 != 3) 
tt3 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

pp(tt3, "L3", "Number of AS molQTLs")

# ggsave("./manuscript_plots/supp/s8.png", width = p_width, height = p_height+0.5)

cp_as1 <- pop2_pip %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05" & L3 != 3) %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L3, cp_as1)) %>%
  mutate(p.value = p.value/2)

cp_as2 <- pop2_cs %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05" & L3 != 3) %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L3, cp_as2)) %>%
  mutate(p.value = p.value/2)

cp_as3 <- pop2_prop %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05" & L3 != 3) %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L3, cp_as3)) %>%
  mutate(p.value = p.value/2)

# L2
pp_l2 <- function(tt1, chvar, xlab, colors_type = method_colors){
  envar = sym(chvar)
  
  p1tmp <- filter(tt1, group == "PIP")
  
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
  
  p2tmp <- filter(tt1, group == "CS")
  
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
  
  p3tmp <- filter(tt1, group == "prec")
  
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

ddPROP2 <- pop2_prop %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L1) %>%
  group_by(sim, name, N, L1, L2, L3, h2g, rho, group) %>%
  summarize(mvalue = mean(value),
    se = 1.96 * (sd(value) / sqrt(n()))) %>%
  mutate(group = "prec")


ttdd <- bind_rows(ddPIP2, ddCS2, ddPROP2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie"))) 

tt4 <- ttdd %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05")

tt4 %>%
  group_by(L1, L2, h2g, rho, L3, N) %>%
  summarize(n = n())

pp_l2(tt4, "L2", "Number of Inferred molQTLs")

# ggsave("./manuscript_plots/supp/s9.png", width = p_width, height = p_height+0.5)

cp_l1 <- pop2_pip %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  select(sushie, indep, meta, susie, sim, locus, N, L1, L2, L3, h2g, rho) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L2, cp_l1)) %>%
  mutate(p.value = p.value/2)

cp_l2 <- pop2_cs %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  select(sushie, indep, meta, susie, sim, CSIndex, locus, N, L1, L2, L3, h2g, rho) %>%
  filter(sushie != 0 & indep != 0 & meta != 0 & susie != 0) %>%
  pivot_longer(cols = sushie:susie) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L2, cp_l2)) %>%
  mutate(p.value = p.value/2)

cp_l3 <- pop2_prop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  pivot_longer(cols = c(prec)) %>%
  rename(group = name,
    name = method) %>%
  mutate(value = value/L2) %>%
  mutate(name = factor(name,
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE"),
    levels = c("sushie", "indep", "meta", "susie")))

tidy(lm(value ~ name + L2, cp_l3)) %>%
  mutate(p.value = p.value/2)

# r2 diff
tmp_r2 <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_pred3_r2.tsv.gz")

df_r2 <- tmp_r2 %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  # filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  # filter(name == "ancestry2_weight2") %>%
  filter(!is.na(value)) %>%
  group_by(sim, N, method) %>%
  summarize(adj_r2 = mean(value),
    adj_se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(method = factor(method,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    N = case_when(
      N == "200:200" ~ 200,
      N == "400:400" ~ 400,
      N == "600:600" ~ 600,
      N == "800:800" ~ 800))

r2_ss <- ggplot(df_r2,
  aes(x = factor(N), y = adj_r2, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = adj_r2 - adj_se, ymax = adj_r2 + adj_se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("Training Sample Size") +
  theme_sim() 


df_r2 <- tmp_r2 %>%
  filter(N == "400:400" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  # filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  # filter(name == "ancestry2_weight2") %>%
  filter(!is.na(value)) %>%
  group_by(sim, h2g, method) %>%
  summarize(adj_r2 = mean(value),
    adj_se = 1.96 * sd(value) / sqrt(n()),
    n = n()) %>%
  mutate(method = factor(method,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2))

r2_h2g <- ggplot(df_r2,
  aes(x = factor(h2g), y = adj_r2, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = adj_r2 - adj_se, ymax = adj_r2 + adj_se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("Predicted r-sqaured") +
  xlab("cis-SNP Heritability") +
  theme_sim() 

cp_r2 <- tmp_r2 %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  # filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  filter(!is.na(value)) %>%
  mutate(method = factor(method,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep", "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")))

tidy(lm(value ~ method + name + N, cp_r2)) %>%
  mutate(p.value = p.value/2)


cp_r2 <- tmp_r2 %>%
  filter(N == "400:400" &
      ngwas == 2e5 & h2ge == 0.3/2000) %>%
  # filter(method %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  pivot_longer(cols=c(ancestry1_weight1, ancestry2_weight2)) %>%
  filter(!is.na(value)) %>%
  mutate(method = factor(method,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")))

tidy(lm(value ~ method + name + h2g, cp_r2)) %>%
  mutate(p.value = p.value/2)


tmp_twas <- read_tsv("~/Documents/github/data/sushie_results/sim/sim_pred3_twas.tsv.gz")

unique(tmp_twas$h2ge)

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
  # filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP"))) %>%
  rename(method = name) %>%
  pivot_longer(cols = c(adj_power, adj_se)) %>%
  separate(name, sep = "_", into = c("type", "name")) %>%
  pivot_wider(values_from = value, names_from = "name")


twas_ss <- ggplot(df_twas,
  aes(x = factor(ngwas), y = power, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = power - se, ymax = power + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("GWAS Sample Size") +
  theme_sim()


df_twas <- tmp_twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & ngwas == 200000) %>%
  pivot_longer(cols=c(sushie, indep, meta, susie, enet, lasso, ridge)) %>%
  mutate(pval = 2 * pnorm(abs(value), lower.tail = FALSE)) %>%
  group_by(1:n()) %>%
  mutate(adjpval = min(1, pval*1000)) %>%
  # filter(ancestry == 2) %>%
  group_by(sim, h2ge, ngwas, name) %>%
  summarize(adj_power = mean(adjpval < 0.05),
    adj_se = 1.96 * sqrt(adj_power * (1 - adj_power) / n())) %>%
  # filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    h2ge = factor(case_when(
      h2ge == unique(tmp_twas$h2ge)[3] ~ "6e-5",
      h2ge == 0.00015 ~ "1.5e-4",
      h2ge == 0.00030 ~ "3e-4",
      h2ge == 0.0006 ~ "6e-4",
    ), levels = c("6e-5",  "1.5e-4", "3e-4", "6e-4"))) %>%
  rename(method = name) %>%
  pivot_longer(cols = c(adj_power, adj_se)) %>%
  separate(name, sep = "_", into = c("type", "name")) %>%
  pivot_wider(values_from = value, names_from = "name")


twas_h2ge <- ggplot(df_twas,
  aes(x = factor(h2ge), y = power, color = method)) +
  geom_point(size=point_size, position=position_dodge(width=0.5)) +
  scale_color_manual(values = twas_colors) +
  geom_errorbar(aes(ymin = power - se, ymax = power + se),
    width = 0.1, position=position_dodge(width=0.5)) +
  scale_y_continuous( labels=scaleFUN) +
  ylab("TWAS Power") +
  xlab("Prop. of trait heritability mediated by expression") +
  theme_sim()

ggarrange(r2_ss, r2_h2g, twas_ss, twas_h2ge,
  common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C", "D"), font.label = list(size = 8))

# ggsave("./manuscript_plots/supp/s10.png", width = p_width, height = p_height+2)

cp_twas <- tmp_twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & h2ge == 0.3/2000) %>%
  pivot_longer(cols=c(sushie, indep, meta, susie, enet, lasso, ridge)) %>%
  mutate(value = value**2) %>%
  # filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
      labels = c("100k", "200k", "300k")))

tidy(lm(value ~ name+ ngwas +factor(ancestry), cp_twas)) %>%
  mutate(p.value = p.value/2)

cp_twas <- tmp_twas %>%
  filter(N == "400:400" & h2g == "0.05:0.05" & ngwas == 200000) %>%
  pivot_longer(cols=c(sushie, indep, meta, susie, enet, lasso, ridge)) %>%
  mutate(value = value**2) %>%
  # filter(name %in% c("sushie", "susie", "lasso", "enet", "ridge")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "indep", "meta", "susie", "lasso", "enet", "ridge"),
    labels = c("SuShiE", "SuShiE-Indep",  "Meta-SuSiE", "SuSiE", "LASSO", "Elastic Net", "gBLUP")),
    ngwas = factor(ngwas,  levels = c(100000, 200000, 300000),
      labels = c("100k", "200k", "300k")))

tidy(lm(value ~ name+ h2ge +factor(ancestry), cp_twas)) %>%
  mutate(p.value = p.value/2)
  

