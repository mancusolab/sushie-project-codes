library(tidyverse)

df_genoa <- read_tsv("/scratch1/zeyunlu/sushie/genoa_sushie_gene_list_noMHC.tsv", col_names = FALSE)
df_geuv <- read_tsv("/scratch1/zeyunlu/sushie/geuvadis_gene_list_noMHC.tsv", col_names = FALSE)

df_same <- df_genoa %>%
  select(ID = X2) %>%
  inner_join(df_geuv %>%
      select(ID=X8),
    by = "ID")

write_tsv(df_same,
  "/scratch1/zeyunlu/sushie/geuvadis_overlap_gene_list_noMHC.tsv", col_names = FALSE)

df_proteins <- read_tsv("/scratch1/zeyunlu/sushie/mesa_proteins_gene_list_noMHC.tsv", col_names = FALSE)
df_interval <- read_tsv("/scratch1/zeyunlu/sushie/interval_gene_list_noMHC.tsv", col_names = FALSE)

df_pro <- df_interval %>%
  select(ID=X8, INTERVAL=X13) %>%
  inner_join(df_proteins %>%
      select(ID=X2, MESA=X15),
    by = "ID", multiple = "all")

write_tsv(df_pro,
  "/scratch1/zeyunlu/sushie/interval_overlap_gene_list_noMHC.tsv", col_names = FALSE)

df_rnaseq <- read_tsv("/scratch1/zeyunlu/sushie/mesa_rnaseq_gene_list_noMHC.tsv", col_names = FALSE)

df_v5 <- read_tsv("/scratch1/zeyunlu/sushie/mesa_rnaseq_v5_gene_list_noMHC.tsv", col_names = FALSE)

valid_v5 <- df_v5 %>%
  select(ID = X12) %>%
  inner_join(df_rnaseq %>%
      select(ID=X12))

write_tsv(valid_v5,
  "/scratch1/zeyunlu/sushie/v5_overlap_gene_list_noMHC.tsv", col_names = FALSE)

# df_tcell <- read_tsv("/scratch1/zeyunlu/sushie/mesa_tcell_gene_list.tsv", col_names = FALSE)
# df_mono <- read_tsv("/scratch1/zeyunlu/sushie/mesa_mono_gene_list.tsv", col_names = FALSE)
# 
# valid_mono <- df_mono %>%
#   select(ID = X11) %>%
#   inner_join(df_rnaseq %>%
#       select(ID=X11))
# 
# valid_tcell <- df_tcell %>%
#   select(ID = X11) %>%
#   inner_join(df_rnaseq %>%
#       select(ID=X11))
# 
# write_tsv(valid_mono,
#   "/scratch1/zeyunlu/sushie/mono_overlap_gene_list.tsv", col_names = FALSE)
# 
# write_tsv(valid_tcell,
#   "/scratch1/zeyunlu/sushie/tcell_overlap_gene_list.tsv", col_names = FALSE)




