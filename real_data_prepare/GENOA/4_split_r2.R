args <- commandArgs(TRUE)
idx <- args[1]

library(tidyverse)

setwd("/project/nmancuso_8/zeyunlu/result")

# EA AFF  1384
# EA ILL 123

tot <- 1384 + 123

fname <- paste0("ea_", idx, "_r2.tsv")
dd <- read_tsv(fname, col_names = FALSE) %>%
  separate(X5, into = c("r1", "r2"), sep =",") %>%
  mutate(r1 = as.numeric(r1),
    r2 = as.numeric(r2),
    rr = ifelse(is.na(r1) & is.na(r2), NA,
      ifelse(is.na(r1) & !is.na(r2), r2,
        ifelse(is.na(r2) & !is.na(r1), r1,
          (1384/tot)*r1 + (123/tot)*r2)))) %>%
  filter(!is.na(rr)) %>%
  select(X1, X2, X3, X4, rr)
write_tsv(dd, paste0("ea_", idx, "_newr2.tsv"), col_names = FALSE)


# AA AFF 1263
# AA ILL 269

tot <- 1263 + 269

fname <- paste0("aa_", idx, "_r2.tsv")
dd <- read_tsv(fname, col_names = FALSE) %>%
  separate(X5, into = c("r1", "r2"), sep =",") %>%
  mutate(r1 = as.numeric(r1),
    r2 = as.numeric(r2),
    rr = ifelse(is.na(r1) & is.na(r2), NA,
      ifelse(is.na(r1) & !is.na(r2), r2,
        ifelse(is.na(r2) & !is.na(r1), r1,
          (1263/tot)*r1 + (269/tot)*r2)))) %>%
  filter(!is.na(rr)) %>%
  select(X1, X2, X3, X4, rr)
write_tsv(dd, paste0("aa_", idx, "_newr2.tsv"), col_names = FALSE)

