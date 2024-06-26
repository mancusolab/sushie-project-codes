# this script is to create main figure 3
library(tidyverse)
library(glue)

library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


# case study URGCP
target_gene <- "ENSG00000106608_URGCP"
target_snp <- "rs2528382"

df_eqtl_eur <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.EUR.qtlscan.PHENO1.glm.linear"))

df_eqtl_afr <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.AFR.qtlscan.PHENO1.glm.linear"))

df_eqtl_his <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.HIS.qtlscan.PHENO1.glm.linear"))

df_eur <- df_eqtl_eur %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target)

df_afr <- df_eqtl_afr %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target)

df_his <- df_eqtl_his %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target)


tmp_weight <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.normal.sushie.weights.tsv"))

df_weight <- tmp_weight %>%
  mutate(chrom = paste0("chr", chrom),
    snp = ifelse(snp %in% target_snp, "", snp),
    sushie_pip_all= sushie_pip_all*0.9) %>%
  dplyr::select(chrom, pos, snp, p = sushie_pip_all, cs = sushie_cs_index)

target_chrom <- "chr7"
target_tss <- 43375894
target_tes <- 43926726 + 500000
snp_pos <- 43926148

df_h3k27ac <- readBigwig(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/URGCP_h3k27ac.bigWig"),
  chrom = target_chrom, chromstart = target_tss, chromend = target_tes)

atac_seq <- read_tsv(
  glue("~/Documents/github/data/sushie_results/case/snATAC-seq-frozen-peaks_grch38.bed"), col_names = FALSE) %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == "chr7" & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

chipseq_params <- pgParams(chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, 
  assembly = "hg38", x = 3.5, width = 1.5, default.units = "inches")

h3k27ac_range <- pgParams(range = c(min(df_h3k27ac$score),
  max(df_h3k27ac$score)), assembly = "hg38")

sn_seq <- read_tsv(glue("~/Documents/github/data/sushie_results/case/chiou_et_al_grch38.bed.gz"), col_names = FALSE) %>%
  separate(X4, sep=",", into=paste0("col", 1:40)) %>%
  pivot_longer(cols = contains("col")) %>%
  filter(value %in% c("cytotoxic-NK", "naive-T", "naive-B",
    "classical-monocyte"))

nk_atac <- sn_seq %>%
  filter(value == "cytotoxic-NK") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == "chr7" & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

t_atac <- sn_seq %>%
  filter(value == "naive-T") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == "chr7" & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

b_atac <- sn_seq %>%
  filter(value == "naive-B") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == "chr7" & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

mono_atac <- sn_seq %>%
  filter(value == "classical-monocyte") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == "chr7" & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))


png(filename="./plots/p4.png", width = 7, height = 5.7, units="in", res=300)

pageCreate(width = 7, height = 5.7, default.units = "inches")

base_height <- 0.5
v_space <- 0.1
super_base <- 0.2
label_font_size <- 8
axis_font_size <- 7
legend_font_size <- 11
figure_width <- 6.2
annot_font_size <- 6

plotText(label = "A", x = 0.1, y = 0.1,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p1 <- plotManhattan(
  data = df_eur, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill = "grey", assembly = "hg38", x = 0.5, y = super_base, width = figure_width,
  height = base_height)

annoYaxis(plot = p1, at = c(0, 1, 2, 3), axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "EUR", x = 6.8, y = super_base + base_height / 2,
  rot = 270, fontsize = label_font_size, fontface = "bold", default.units = "inches")

p2 <- plotManhattan(
  data = df_afr, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"), fill="grey",
  assembly = "hg38", x = 0.5, y = super_base + v_space + base_height,
  width = figure_width, height = base_height)

annoYaxis(plot = p2, at = c(0, 1, 2, 3, 4), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "AFR", x = 6.8,
  y = super_base + v_space + base_height + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

p3 <- plotManhattan(
  data = df_his, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"), fill="grey",
  assembly = "hg38", x = 0.5,
  y = super_base + v_space*2 + base_height*2,
  width = figure_width, height = base_height
)

annoYaxis(plot = p3, at = c(0, 1, 2, 3, 4), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "HIS", x = 6.8,
  y = super_base + v_space*2 + base_height*2 + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

plotText(label = "-log10(P)", x = 0.2,
  y = super_base + (v_space*2 + base_height*3) / 2, rot = 90,
  fontsize = annot_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p4 <- plotManhattan(data = df_weight, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, trans = "",
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill="grey", assembly = "hg38", x = 0.5,
  y = super_base + v_space*3 + base_height*3,
  width = figure_width, height = base_height*0.75)

annoYaxis(plot = p4, at = c(0, 1), axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "PIP", x = 0.2,
  y = super_base + v_space*3 + base_height*3 + (base_height*0.75) / 2,
  rot = 90, fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p5 <- plotGenes(assembly = "hg38", chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes, geneHighlights = data.frame("gene" = c("URGCP"),
    "color" = c("orange")),
  geneBackground = "grey", fontsize = axis_font_size, x = 0.5,
  y = super_base + v_space*4 + base_height*3.75,
  width=figure_width, height = base_height)

plotText(label = "B", x = 0.1, y = super_base + v_space*4 + base_height*4.75,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p6 <- plotSignal(data = df_h3k27ac,
  params = c(chipseq_params, h3k27ac_range),
  fill = "#253494", linecolor = "#253494",
  x = 0.5, y = super_base + v_space*4 + base_height*4.75 + 0.2,
  height = base_height*0.75, width = figure_width)

annoYaxis(plot = p6, at = c(0, 40), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "H3K27ac", x = 0.2,
  y = super_base + v_space*4 + base_height*4.75 + 0.2 + (base_height*0.75) / 2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", 
  default.units = "inches")

p7 <- plotRanges(
  data = atac_seq, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = super_base + v_space*5 + base_height*5.5 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "PBMC", x = 0.2,
  y = super_base + v_space*5 + base_height*5.5 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p8 <- plotRanges(
  data = t_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = super_base + v_space*6 + base_height*6 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Naive T", x = 0.2,
  y = super_base + v_space*6 + base_height*6 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p9 <- plotRanges(
  data = b_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = super_base + v_space*7 + base_height*6.5 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Naive B", x = 0.2,
  y = super_base + v_space*7 + base_height*6.5 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p10 <- plotRanges(
  data = nk_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = super_base + v_space*8 + base_height*7 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "cNK", x = 0.2,
  y = super_base + v_space*8 + base_height*7 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p11 <- plotRanges(
  data = mono_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = super_base + v_space*9 + base_height*7.5 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Monocyte", x = 0.2,
  y = super_base + v_space*9 + base_height*7.5 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

annoGenomeLabel(plot= p8, x = 0.5,
  y = super_base + v_space*9 + base_height*8.25 + 0.2,
  fontsize = label_font_size, scale = "Mb", default.units = "inches")

pageGuideHide()

dev.off()

