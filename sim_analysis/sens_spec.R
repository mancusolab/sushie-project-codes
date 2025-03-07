library(tidyverse)
library(ggpubr)
library(broom)
library(RColorBrewer)

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

df2 <- crossing(locus=1:500, name = c("SuShiE", "SuShiE-Indep", "SuSiE",
  "Meta-SuSiE",
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

# average sens
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

# average specificity 
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

# average FDR
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

