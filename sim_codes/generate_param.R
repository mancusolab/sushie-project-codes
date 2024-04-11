library(tidyverse)

num_gene <- 500

gene_list <- read_tsv("sim_gene_list.tsv", col_names = FALSE) %>%
  separate(col = X1, into=c("pop", "gene", "name", "suffix"), sep = "_") %>%
  mutate(LOCUS = as.numeric(gsub("gene", "", gene))) %>%
  select(gene = name, LOCUS) %>%
  filter(LOCUS <= num_gene)


# 2 pop
# N
N <- c("200:200", "400:400", "600:600", "800:800")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt1 <- crossing(N, L1, L2, L3, h2g, rho)

# L1
N <- "400:400"
L1 <- c(1, 2, 3)
L2 <- c(1, 2, 3)
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt2 <- crossing(N, L1, L2, L3, h2g, rho) %>%
  filter(L1 == L2)

# h2g
N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")
rho <- "0.8"
tt3 <- crossing(N, L1, L2, L3, h2g, rho)

# rho
N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- c("0.01", "0.4", "0.8", "0.99")
tt4 <- crossing(N, L1, L2, L3, h2g, rho)




# model assumption unmet
# different N
N <- c("400:200", "400:400", "400:600", "400:800")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt5 <- crossing(N, L1, L2, L3, h2g, rho)

# different h2g
N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")
rho <- "0.8"
tt6 <- crossing(N, L1, L2, L3, h2g, rho)

# L2 1
N <- "400:400"
L1 <- c(1, 2, 3)
L2 <- 5
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt7 <- crossing(N, L1, L2, L3, h2g, rho)

# L2 2
N <- "400:400"
L1 <- 2
L2 <- c(2, 5, 10)
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt8 <- crossing(N, L1, L2, L3, h2g, rho)

# L3
N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- c(0, 1, 2)
h2g <- "0.05:0.05"
rho <- "0.8"
tt9 <- crossing(N, L1, L2, L3, h2g, rho)

# rho performance
N <- c("400:400", "1200:1200", "2400:2400")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
tt10 <- crossing(N, L1, L2, L3, h2g, rho)

tt <- bind_rows(tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, tt10) %>%
  distinct() %>%
  rownames_to_column("ROW") 

dd1 <- tt %>%
  left_join(crossing(ROW = tt$ROW, LOCUS = 1:500),
    by = "ROW") %>%
  mutate(seed = row_number() + 100) %>%
  left_join(gene_list, by = "LOCUS")

write_tsv(dd1, file="./param/param_2pop.tsv", col_names = FALSE, escape = "none")


# 3 pop
N <- c("450", "900", "1350")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05"
rho <- "0.8"
tt1 <-crossing(N, L1, L2, L3, h2g, rho)


tt <- bind_rows(tt1) %>%
  distinct() %>%
  rownames_to_column("ROW") 

dd2 <- tt %>%
  left_join(crossing(ROW = tt$ROW, LOCUS = 1:500),
    by = "ROW") %>%
  mutate(seed = row_number() + 100) %>%
  left_join(gene_list, by = "LOCUS")

write_tsv(dd2, file="./param/param_3pop.tsv", col_names = FALSE, escape = "none")


# no shared
N <- c("400:200", "400:400", "400:600", "400:800")
L1 <- 0
L2 <- 2
L3 <- 2
h2g <- "0.05:0.05"
rho <- "0.8"
tt1 <- crossing(N, L1, L2, L3, h2g, rho)

N <- "400:400"
L1 <- 0
L2 <- 2
L3 <- 2
h2g <- c("0.05:0.01", "0.05:0.05", "0.05:0.1", "0.05:0.2")
rho <- "0.8"
tt2 <- crossing(N, L1, L2, L3, h2g, rho)

tt <- bind_rows(tt1, tt2) %>%
  distinct() %>%
  rownames_to_column("ROW") 

dd3 <- tt %>%
  left_join(crossing(ROW = tt$ROW, LOCUS = 1:500),
    by = "ROW") %>%
  mutate(seed = row_number() + 100)  %>%
  left_join(gene_list, by = "LOCUS")

write_tsv(dd3, file="./param/param_noshared.tsv", col_names = FALSE, escape = "none")


# pred

N <- c("200:200", "400:400", "600:600")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
ngwas <- 2e5
h2ge <- 0.3/2000
tt1 <- crossing(N, L1, L2, L3, h2g, rho, ngwas, h2ge)


N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- c("0.01:0.01", "0.05:0.05", "0.1:0.1", "0.2:0.2")
rho <- "0.8"
ngwas <- 2e5
h2ge <- 0.3/2000
tt2 <- crossing(N, L1, L2, L3, h2g, rho, ngwas, h2ge)

N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
ngwas <- c(1e5, 2e5, 3e5)
h2ge <- 0.3/2000
tt3 <- crossing(N, L1, L2, L3, h2g, rho, ngwas, h2ge)

N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- "0.8"
ngwas <- 2e5
h2ge <- c(0.3/5000, 0.3/2000, 0.3/1000, 0.3/500)
tt4 <- crossing(N, L1, L2, L3, h2g, rho, ngwas, h2ge)

tt <- bind_rows(tt1, tt2, tt3, tt4) %>%
  distinct() %>%
  rownames_to_column("ROW") 

dd4 <- tt %>%
  left_join(crossing(ROW = tt$ROW, LOCUS = 1:500),
    by = "ROW") %>%
  mutate(seed = row_number() + 100)  %>%
  left_join(gene_list, by = "LOCUS")

write_tsv(dd4, file="./param/param_pred.tsv", col_names = FALSE, escape = "none")


# rho
# rho performance
N <- c("400:400", "800:800")
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- "0.05:0.05"
rho <- c("0.01", "0.4", "0.8", "0.99")
tt1 <- crossing(N, L1, L2, L3, h2g, rho)

N <- "400:400"
L1 <- c(2, 3)
L2 <- c(2, 3)
L3 <- 0
h2g <- "0.05:0.05"
rho <- c("0.01", "0.4", "0.8", "0.99")
tt2 <- crossing(N, L1, L2, L3, h2g, rho) %>%
  filter(L1==L2)

N <- "400:400"
L1 <- 2
L2 <- 2
L3 <- 0
h2g <- c("0.05:0.05", "0.2:0.2")
rho <- c("0.01", "0.4", "0.8", "0.99")
tt3 <- crossing(N, L1, L2, L3, h2g, rho)

tt <- bind_rows(tt1, tt2, tt3) %>%
  distinct() %>%
  rownames_to_column("ROW") 

dd5 <- tt %>%
  left_join(crossing(ROW = tt$ROW, LOCUS = 1:500),
    by = "ROW") %>%
  mutate(seed = row_number() + 100) %>%
  left_join(gene_list, by = "LOCUS")

write_tsv(dd5, file="./param/param_rho.tsv", col_names = FALSE, escape = "none")


