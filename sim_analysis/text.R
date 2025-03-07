library(tidyverse)
library(broom)
library(glue)

sim_data_path <- "~/Downloads/sushie_sim_data_results/"

other_methods <- c("SuShiE-Indep",
  "Meta-SuSiE", "SuSiE", "SuSiEx", "MESuSiE", 
  "XMAP", "XMAP-IND")

# 3 pop
load("./data/df_3pop.RData")
comp_methods_3pop <- function(df, method, type_name) {
  tmp_cp <- df %>%
    filter(type %in% type_name) %>%
    filter(name %in% c("Three-Ancestry", method)) %>%
    mutate(name = factor(name, levels = c("Three-Ancestry", method))) %>%
    select(-type, -sim, -locus)
  
  tmp_res <- tidy(lm(value ~ name + ., tmp_cp)) %>%
    filter(grepl("name", term)) %>%
    mutate(term = gsub("name", "", term)) %>%
    select(term, estimate, se = std.error, p.value)
  return(tmp_res)
}

cp_res <- tibble()
for (other_method in c("One-Ancestry", "Two-Ancestry")) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods_3pop(df_3pop, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# finish inference
sushie_sens_pop2 <- read_tsv(glue("{sim_data_path}/sushie_2pop_sens.tsv.gz"))
mesusie_sens_pop2 <- read_tsv(glue("{sim_data_path}/mesusie_in_2pop_sens.tsv.gz"))
xmap_sens_pop2 <- read_tsv(glue("{sim_data_path}/xmap_in_2pop_sens.tsv.gz"))
xmap_ind_sens_pop2 <- read_tsv(glue("{sim_data_path}/xmap_ind_2pop_sens.tsv.gz"))
susiex_sens_pop2 <- read_tsv(glue("{sim_data_path}/susiex_in_2pop_sens.tsv.gz"))

df_sens <- sushie_sens_pop2 %>%
  pivot_longer(cols = c(sushie:susie)) %>%
  select(name, value, sim, locus) %>%
  bind_rows(mesusie_sens_pop2 %>%
      mutate(name = "mesusie") %>%
      rename(value = mesusie),
    xmap_sens_pop2 %>%
      mutate(name = "xmap") %>%
      rename(value = xmap),
    xmap_ind_sens_pop2 %>%
      mutate(name = "xmap_ind") %>%
      rename(value = xmap),
    susiex_sens_pop2 %>%
      rename(value = sens,
        name = method)
  ) %>%
  group_by(name) %>%
  summarize(mvalue = mean(value),
    se = sd(value)/sqrt(n())) %>%
  filter(name %in% c("sushie", "mesusie", "xmap", "xmap_ind", "susiex")) %>%
  mutate(name = factor(name,
    levels = c("sushie", "mesusie", "xmap", "xmap_ind", "susiex"),
    labels = c("SuShiE", "MESuSiE", "XMAP", "XMAP-IND", "SuSiEx"))) %>%
  arrange(name)

(0.255 + 0.261)/2
# write_tsv(df_sens, "./manuscript_plots/s2.tsv")

# 2 pop general data
load("./data/df_2pop.RData")

comp_methods <- function(df, method, type_name) {
  tmp_cp <- df %>%
    filter(type %in% type_name) %>%
    filter(name %in% c("SuShiE", method)) %>%
    mutate(name = factor(name, levels = c("SuShiE", method))) %>%
    select(-type, -sim, -locus)
  
  tmp_res <- tidy(lm(value ~ name + ., tmp_cp)) %>%
    filter(grepl("name", term)) %>%
    mutate(term = gsub("name", "", term)) %>%
    select(term, estimate, se = std.error, p.value)
  return(tmp_res)
}

# across all normal simulations, no misspecification
df_tmp <- df_2pop %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == L2 & L3 == 0 &
      h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800),
    h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2)) %>%
  select(-L3, -L2) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  filter(grepl("XMAP", term)) %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

xmap_only_case <- df_tmp %>%
  filter(type == "CS") %>%
  filter(name %in% c("SuShiE", "XMAP")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(!is.na(XMAP)) %>%
  pivot_longer(cols = c(SuShiE:XMAP))

xmap_ind_only_case <- df_tmp %>%
  filter(type == "CS") %>%
  filter(name %in% c("SuShiE", "XMAP-IND")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(!is.na(`XMAP-IND`)) %>%
  pivot_longer(cols = c(SuShiE:`XMAP-IND`))

df_tmp %>%
  filter(type == "CS") %>%
  filter(name %in% "SuShiE") %>%
  filter(!is.na(value)) %>%
  distinct(sim, locus, CSIndex) %>%
  left_join(df_tmp %>%
      filter(type == "CS") %>%
      filter(name %in% c("SuShiE", "XMAP")) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      filter(!is.na(XMAP)) %>%
      pivot_longer(cols = c(SuShiE:XMAP)) %>%
      distinct(sim, locus, CSIndex) %>%
      mutate(XMAP = 1)) %>%
  left_join(
    df_tmp %>%
      filter(type == "CS") %>%
      filter(name %in% c("SuShiE", "XMAP-IND")) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      filter(!is.na(`XMAP-IND`)) %>%
      pivot_longer(cols = c(SuShiE:`XMAP-IND`)) %>%
      distinct(sim, locus, CSIndex) %>%
      mutate(`XMAP-IND` = 1)) %>%
  ungroup() %>%
  summarize(prop = sum(!(is.na(XMAP) & is.na(`XMAP-IND`)))/n())

bind_rows(comp_methods(xmap_only_case, "XMAP", "CS"),
  comp_methods(xmap_ind_only_case, "XMAP-IND", "CS")) %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

# main figure 2
df_tmp <- df_2pop %>%
  filter(N %in% c("200:200", "400:400", "600:600", "800:800") &
      L1 == L2 & L3 == 0 &
      h2g %in%  "0.05:0.05" & 
      rho %in%  0.8) %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800)) %>%
  select(-L3, -L2, -h2g, -L1, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# Figure S2-S4
df_tmp <- df_2pop %>%
  filter(N %in% "400:400" &
      L1 == L2 & L3 == 0 &
      h2g %in% c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")) %>%
  mutate(h2g = case_when(
    h2g == "0.01:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.1:0.1" ~ 0.1,
    h2g == "0.2:0.2" ~ 0.2)) %>%
  select(-L3, -L2, -N, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# differet N
df_tmp <- df_2pop %>%
  filter(L1 == 2 & L2 ==2 & L3 == 0 & h2g == "0.05:0.05" & rho == 0.8) %>%
  filter(N %in% c("400:200", "400:400", "400:600", "400:800")) %>%
  mutate(N = case_when(
    N == "400:200" ~ 200,
    N == "400:400" ~ 400,
    N == "400:600" ~ 600,
    N == "400:800" ~ 800)) %>%
  select(-L3, -L2, -L1, -h2g, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS")

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

# different h2g

df_tmp <- df_2pop %>%
  filter(N == "400:400" & L1 == 2 & L2 ==2 & L3 == 0 & rho == 0.8) %>%
  filter(h2g %in% c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")) %>%
  mutate(h2g = case_when(
    h2g == "0.05:0.01" ~ 0.01,
    h2g == "0.05:0.05" ~ 0.05,
    h2g == "0.05:0.1" ~ 0.1,
    h2g == "0.05:0.2" ~ 0.2))  %>%
  select(-L3, -L2, -L1, -N, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# additional AS QTL
df_tmp <- df_2pop %>%
  filter(N == "400:400" & L2 == 2 &  L1 ==2  &
      rho == 0.8 & h2g == "0.05:0.05") %>%
  select(-h2g, -L2, -L1, -N, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


# additional L, comparison
df_tmp <- df_2pop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  filter(!is.na(value)) %>%
  select(-h2g, -L1, -L3, -N, -rho) 

cp_res <- tibble()
for (other_method in other_methods) {
  for (type_name in c("PIP", "CS", "Calibration")) {
    cp_res <- cp_res %>%
      bind_rows(
        comp_methods(df_tmp, other_method, type_name) %>%
          mutate(type = type_name) %>%
          mutate(mval = mean(estimate))
      )
  }
}

cp_res %>%
  filter(type == "PIP") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

cp_res %>%
  filter(type == "CS") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    weighted_se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/weighted_se), lower.tail = FALSE))

cp_res %>%
  filter(type == "Calibration") %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

# additional L, L's effect 
tt5 <- df_2pop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  filter(!is.na(value)) %>%
  filter(name %in% "SuShiE") %>%
  filter(type %in% "PIP")

tidy(lm(value ~ L2 + CSIndex, data = tt5))

tt6 <- df_2pop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  filter(!is.na(value)) %>%
  filter(name %in% "SuShiE") %>%
  filter(type %in% "CS")

tidy(lm(value ~ L2 + CSIndex, data = tt6))

tt7 <- df_2pop %>%
  filter(N == "400:400" &  L1 ==2 & L3 == 0 & rho == 0.8 & h2g == "0.05:0.05") %>%
  filter(!is.na(value)) %>%
  filter(name %in% "SuShiE") %>%
  filter(type %in% "Calibration")

tidy(lm(value ~ L2 + CSIndex, data = tt7))

# r2 and twas
load("./data/df_r2twas.RData")

comp_twas <- function(df, method) {
  tmp_cp <- df %>%
    filter(name %in% c("SuShiE", method)) %>%
    mutate(name = factor(name, levels = c("SuShiE", method))) %>%
    select(-sim, -locus)
  
  tmp_res <- tidy(lm(value ~ name + ., tmp_cp)) %>%
    filter(grepl("name", term)) %>%
    mutate(term = gsub("name", "", term)) %>%
    select(term, estimate, se = std.error, p.value)
  return(tmp_res)
}

# main figure
df_tmp <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000 & type == "R2") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800)) %>%
  select(-ngwas, -h2ge, -h2g, -type)

cp_res <- tibble()
for (other_method in other_methods[!other_methods %in% "SuSiEx"]) {
  cp_res <- cp_res %>%
    bind_rows(
      comp_twas(df_tmp, other_method) %>%
        mutate(mval = mean(estimate))
    )
}

cp_res %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

df_tmp <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      N == "400:400" & h2ge == 0.3/2000 & type == "TWAS") %>%
  mutate(ngwas = factor(ngwas,  levels = c(100000, 200000, 300000))) %>%
  select(-N, -h2ge, -h2g, -type)

cp_res <- tibble()
for (other_method in other_methods[!other_methods %in% "SuSiEx"]) {
  cp_res <- cp_res %>%
    bind_rows(
      comp_twas(df_tmp, other_method) %>%
        mutate(mval = mean(estimate))
    )
}

cp_res %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

# supp figure
df_tmp <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      ngwas == 2e5 & h2ge == 0.3/2000 & type == "R2") %>%
  mutate(N = case_when(
    N == "200:200" ~ 200,
    N == "400:400" ~ 400,
    N == "600:600" ~ 600,
    N == "800:800" ~ 800),
    h2g = case_when(
      h2g == "0.01:0.01" ~ 0.01,
      h2g == "0.05:0.05" ~ 0.05,
      h2g == "0.1:0.1" ~ 0.1,
      h2g == "0.2:0.2" ~ 0.2)) %>%
  select(-ngwas, -h2ge, -type)

cp_res <- tibble()
for (other_method in c(other_methods[!other_methods %in% "SuSiEx"],
  "Elastic Net", "LASSO", "gBLUP")) {
  cp_res <- cp_res %>%
    bind_rows(
      comp_twas(df_tmp, other_method) %>%
        mutate(mval = mean(estimate))
    )
}

cp_res %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))

df_tmp <- df_r2twas %>%
  filter(h2g == "0.05:0.05" &
      N == "400:400" & type == "TWAS") %>%
  select(-N, -h2g, -type)

cp_res <- tibble()
for (other_method in c(other_methods[!other_methods %in% "SuSiEx"],
  "Elastic Net", "LASSO", "gBLUP")) {
  cp_res <- cp_res %>%
    bind_rows(
      comp_twas(df_tmp, other_method) %>%
        mutate(mval = mean(estimate))
    )
}

cp_res %>%
  mutate(weight = 1/(se^2)) %>%
  summarize(weighted_mean = sum(estimate * weight) / sum(weight),
    se = sqrt(1/sum(weight)),
    p.value = 2*pnorm(abs(weighted_mean/se), lower.tail = FALSE))


load("./data/df_sushie_comp.RData")

cor(df_sushie_comp$sushie, df_sushie_comp$sushie_ss)

