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

load("./data/df_2pop.RData")

df_cali <- df_2pop %>%
  filter(type == "Calibration") %>%
  # filter(type == "CS") %>%
  filter(N== "400:400" & L1 == 1 & L2 == 1 & L3 == 0 & rho == "0.8") %>%
  filter(h2g %in% c("0.05:0.05")) %>%
  ungroup() %>%
  select(name, locus, value)

df1 <- crossing(locus=1:500,
  name = c("SuShiE", "SuShiE-Indep", "SuSiE", "Meta-SuSiE",
  "SuSiEx", "MESuSiE")) %>%
  mutate(TR = 1) %>%
  left_join(df_cali,
    by = c("locus", "name"))

df_cs <- df_2pop %>%
  filter(type == "CS") %>%
  filter(h2g %in% c("0:0"))

df2 <- crossing(locus=1:500, name = c("SuShiE", "SuShiE-Indep", "SuSiE", "Meta-SuSiE",
  "SuSiEx", "MESuSiE")) %>%
  mutate(TR = 0) %>%
  left_join(df_cs %>%
      mutate(value = 1) %>%
      dplyr::select(name, locus, value),
    by = c("locus", "name")) %>%
  mutate(value = ifelse(is.na(value), 0, 1))

df_total <- bind_rows(df1, df2) %>%
  mutate(value = ifelse(value == 0, "Negative", "Positive"),
    TR = ifelse(TR == 1, "True", "False"))

df_total %>%
  group_by(name) %>%
  summarize(TP = sum(TR == "True" & value == "Positive"),
    FP = sum(TR == "False" & value == "Positive"),
    TN = sum(TR == "False" & value == "Negative"),
    FN = sum(TR == "True" & value == "Negative")) %>%
  mutate(sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    FPR = FP / (FP + TN),
    FDR = FP / (TP + FP))

df_total %>%
  group_by(name) %>%
  summarize(TP = sum(TR == "True" & value == "Positive"),
    FP = sum(TR == "False" & value == "Positive"),
    TN = sum(TR == "False" & value == "Negative"),
    FN = sum(TR == "True" & value == "Negative")) %>%
  mutate(sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    FPR = FP / (FP + TN),
    FDR = FP / (TP + FP)) %>%
  filter(name != "SuShiE") %>%
  summarize(mval = mean(sens))

df_total %>%
  group_by(name) %>%
  summarize(TP = sum(TR == "True" & value == "Positive"),
    FP = sum(TR == "False" & value == "Positive"),
    TN = sum(TR == "False" & value == "Negative"),
    FN = sum(TR == "True" & value == "Negative")) %>%
  mutate(sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    FPR = FP / (FP + TN),
    FDR = FP / (TP + FP))%>%
  filter(name != "SuShiE") %>%
  summarize(mval = mean(spec))

df_total %>%
  group_by(name) %>%
  summarize(TP = sum(TR == "True" & value == "Positive"),
    FP = sum(TR == "False" & value == "Positive"),
    TN = sum(TR == "False" & value == "Negative"),
    FN = sum(TR == "True" & value == "Negative")) %>%
  mutate(sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    FPR = FP / (FP + TN),
    FDR = FP / (TP + FP))%>%
  filter(name != "SuShiE") %>%
  summarize(mval = mean(FPR))

df_total %>%
  group_by(name) %>%
  summarize(TP = sum(TR == "True" & value == "Positive"),
    FP = sum(TR == "False" & value == "Positive"),
    TN = sum(TR == "False" & value == "Negative"),
    FN = sum(TR == "True" & value == "Negative")) %>%
  mutate(sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    FPR = FP / (FP + TN),
    FDR = FP / (TP + FP))%>%
  filter(name != "SuShiE") %>%
  summarize(mval = mean(FDR))

