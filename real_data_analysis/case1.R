# this script is to create main figure 3
library(tidyverse)
library(glue)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(RColorBrewer)

# data folder
case_data_folder <- "~/Downloads/sushie_real_data_results/case_study/"

# case study URGCP
target_gene <- "ENSG00000106608_URGCP"
target_snp <- "rs2528382"

# check if this gene has credible sets in other results
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_normal.sushie_cs.tsv.gz
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_indep.sushie_cs.tsv.gz
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_normal.meta_cs.tsv.gz
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_normal.mega_cs.tsv.gz
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_susiex_cs.tsv.gz
# zgrep "URGCP" ~/Documents/github/data/sushie_results/real2/rnaseq_mesusie_cs.tsv.gz

df_eqtl_eur <- read_tsv(glue("{case_data_folder}/{target_gene}/{target_gene}.EUR.qtlscan.PHENO1.glm.linear"))

df_eqtl_afr <- read_tsv(glue("{case_data_folder}/{target_gene}/{target_gene}.AFR.qtlscan.PHENO1.glm.linear"))

df_eqtl_his <- read_tsv(glue("{case_data_folder}/{target_gene}/{target_gene}.HIS.qtlscan.PHENO1.glm.linear"))

dd <- df_eqtl_eur %>%
  dplyr::select(CHROM = `#CHROM`, POS, ID, REF, ALT,
    EUR_BETA = BETA, EUR_SE = SE, EUR_T_STAT = T_STAT) %>%
  left_join(df_eqtl_afr %>%
      dplyr::select(ID, AFR_BETA = BETA, AFR_SE = SE, AFR_T_STAT = T_STAT)) %>%
  left_join(df_eqtl_his %>%
      dplyr::select(ID, HIS_BETA = BETA, HIS_SE = SE, HIS_T_STAT = T_STAT))

# write_tsv(dd, "./tables/s10.tsv")

target_chrom <- "chr7"
target_tss <- 43375894
target_tes <- 44426410
snp_pos <- 43926148

tmp_weight <- read_tsv(glue(
  "{case_data_folder}/{target_gene}/{target_gene}.normal.sushie.weights.tsv"))

snp_index <- which(tmp_weight$snp == target_snp)

snp_r2_1 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.ans0.{target_snp}.ldr2.tsv"))

snp_r2_2 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.ans1.{target_snp}.ldr2.tsv"))


snp_r2_3 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.ans2.{target_snp}.ldr2.tsv"))

df_eur <- df_eqtl_eur %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target) %>%
  left_join(snp_r2_1) %>%
  mutate(LDgrp = addNA(cut(ldr2, c(0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0)))) %>%
  arrange(LDgrp)

df_afr <- df_eqtl_afr %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target) %>%
  left_join(snp_r2_2) %>%
  mutate(LDgrp = addNA(cut(ldr2, c(0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0)))) %>%
  arrange(LDgrp)

df_his <- df_eqtl_his %>%
  mutate(chrom = paste0("chr", `#CHROM`),
    target = ifelse(ID %in% target_snp, 1, 0),
    ID = ifelse(ID %in% target_snp, "", ID)) %>%
  dplyr::select(chrom, pos = POS, p =P, snp = ID, target) %>%
  left_join(snp_r2_3) %>%
  mutate(LDgrp = addNA(cut(ldr2, c(0, 0.15, 0.30, 0.45, 0.6, 0.75, 1.0)))) %>%
  arrange(LDgrp)

df_weight <- tmp_weight %>%
  mutate(chrom = paste0("chr", chrom),
    snp = ifelse(snp %in% target_snp, "", snp),
    sushie_pip_all= sushie_pip_all*0.9) %>%
  dplyr::select(chrom, pos, snp, p = sushie_pip_all, cs = sushie_cs_index)


df_h3k27ac <- readBigwig(glue(
  "{case_data_folder}/{target_gene}/URGCP_h3k27ac.bigWig"),
  chrom = target_chrom, chromstart = target_tss, chromend = target_tes)

atac_seq <- read_tsv(
  glue("{case_data_folder}/snATAC-seq-frozen-peaks_grch38.bed"), col_names = FALSE) %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

chipseq_params <- pgParams(chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, 
  assembly = "hg38", x = 3.5, width = 1.5, default.units = "inches")

h3k27ac_range <- pgParams(range = c(min(df_h3k27ac$score),
  max(df_h3k27ac$score)), assembly = "hg38")

sn_seq <- read_tsv(glue("{case_data_folder}/chiou_et_al_grch38.bed.gz"), col_names = FALSE) %>%
  separate(X4, sep=",", into=paste0("col", 1:40)) %>%
  pivot_longer(cols = contains("col")) %>%
  filter(value %in% c("cytotoxic-NK", "naive-T", "naive-B",
    "classical-monocyte"))

nk_atac <- sn_seq %>%
  filter(value == "cytotoxic-NK") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

t_atac <- sn_seq %>%
  filter(value == "naive-T") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

b_atac <- sn_seq %>%
  filter(value == "naive-B") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

mono_atac <- sn_seq %>%
  filter(value == "classical-monocyte") %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

df_ld1 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.0.trueLD.tsv")) %>%
  filter(counts != 1)

df_ld2 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.1.trueLD.tsv")) %>%
  filter(counts != 1) 

df_ld3 <- read_tsv(
  glue("{case_data_folder}/{target_gene}/{target_gene}.2.trueLD.tsv")) %>%
  filter(counts != 1)  
max_counts <- max(df_ld1$counts, df_ld2$counts, df_ld3$counts)

df_pairs <- tibble("chrom1" = target_chrom, "start1" = snp_pos, "end1" = snp_pos,
  "chrom2" = target_chrom, "start2" = snp_pos, "end2" = snp_pos)

png(filename="./plots/p4.png", width = 7, height = 7.2, units="in", res=300)

base_height <- 0.5
v_space <- 0.1
super_base <- 0.2
label_font_size <- 8
axis_font_size <- 7
legend_font_size <- 11
figure_width <- 6.2
annot_font_size <- 6

pageCreate(width = 7, height = 7.2, default.units = "inches")

plotText(label = "A", x = 0.1, y = 0.1,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p_ld1 <- plotHicTriangle(df_ld1, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0, max_counts),
  x=0.6, y=super_base , width = 1.8,
  height = base_height*1.8, just = c("left", "top"),
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld1, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "EUR", x = 0.6+1.8/2, y =super_base + base_height*1.8 +0.1,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

p_ld2<- plotHicTriangle(df_ld2, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0,max_counts),
  x=0.6 + 1.8 + 0.25, y=super_base , width = 1.8,
  height = base_height*1.8, just = c("left", "top"),
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld2, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "AFR", x = 0.6 + 1.8 + 0.25 +1.8/2,
  y =super_base + base_height*1.8 +0.1,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

p_ld3 <- plotHicTriangle(df_ld3, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0,max_counts),
  x=0.6 + 1.8*2 + 0.25*2, y=super_base , width = 1.8,
  height = base_height*1.8, just = c("left", "top"),
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld3, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "HIS", x = 0.6 + 1.8*2 + 0.25*2 + 1.8/2,
  y =super_base + base_height*1.8 +0.1,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

annoHeatmapLegend(
  plot = p_ld3,
  x = 0.6 + 1.8*3 + 0.25*2 + 0.1, y = super_base,
  width = 0.2, height = base_height*1.8,
  just = c("left", "top"), default.units = "inches"
)

p1 <- plotManhattan(
  data = df_eur, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes, 
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill=colorby("LDgrp",
    palette = colorRampPalette(c(
      "grey", "#66c2a4",
      "#41ae76", "#238b45", "#006d2c",
      "#00441b", "red"
    ))),
  assembly = "hg38", x = 0.5,
  y = 1.3+super_base, width = figure_width,
  height = base_height)

annoYaxis(plot = p1, at = c(0, 1, 2, 3), axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "EUR", x = 6.8, y = 1.3+super_base + base_height / 2,
  rot = 270, fontsize = label_font_size, fontface = "bold", default.units = "inches")

p2 <- plotManhattan(
  data = df_afr, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill=colorby("LDgrp",
    palette = colorRampPalette(c(
      "grey", "#66c2a4",
      "#41ae76", "#238b45", "#006d2c",
      "#00441b", "red"
    ))), 
  assembly = "hg38", x = 0.5, y = 1.3+super_base + v_space + base_height,
  width = figure_width, height = base_height)

annoYaxis(plot = p2, at = c(0, 1, 2, 3, 4), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "AFR", x = 6.8,
  y = 1.3+super_base + v_space + base_height + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

p3 <- plotManhattan(
  data = df_his, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill=colorby("LDgrp",
    palette = colorRampPalette(c(
      "grey", "#66c2a4",
      "#41ae76", "#238b45", "#006d2c",
      "#00441b", "red"
    ))), 
  assembly = "hg38", x = 0.5,
  y = 1.3+super_base + v_space*2 + base_height*2,
  width = figure_width, height = base_height
)

annoYaxis(plot = p3, at = c(0, 1, 2, 3, 4), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "HIS", x = 6.8,
  y = 1.3+super_base + v_space*2 + base_height*2 + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

plotLegend(
  legend = c(
    target_snp,
    "0-0.15",
    "0.15-0.3",
    "0.3-0.45",
    "0.45-0.6",
    "0.6-0.75",
    "0.75-1"
  ),
  fill = c("red", "grey", "#66c2a4", "#41ae76", "#238b45",
    "#006d2c", "#00441b"),
  cex = 0.9,
  pch = c(18, 19, 19, 19,19,19,19), border = FALSE, x = 2, y = 3.3,
  width = 1.5, height = 0.6, fontface = "bold",just = "center",
  default.units = "inches", orientation = "h", fontsize = label_font_size
)

plotText(label = "-log10(P)", x = 0.2,
  y = 1.3+super_base + (v_space*2 + base_height*3) / 2, rot = 90,
  fontsize = annot_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p4 <- plotManhattan(data = df_weight, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, trans = "",
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill="grey", assembly = "hg38", x = 0.5,
  y = 1.3+super_base + v_space*3 + base_height*3+0.2,
  width = figure_width, height = base_height*0.75)

annoYaxis(plot = p4, at = c(0, 1), axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "PIP", x = 0.2,
  y = 1.3+super_base + v_space*3 + base_height*3 + (base_height*0.75) / 2+0.2,
  rot = 90, fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p5 <- plotGenes(assembly = "hg38", chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes, geneHighlights = data.frame("gene" = c("URGCP"),
    "color" = c("orange")),
  geneBackground = "grey", fontsize = axis_font_size, x = 0.5,
  y = 1.3+super_base + v_space*4 + base_height*3.75 +0.2,
  width=figure_width, height = base_height)

plotText(label = "B", x = 0.1, y = 1.3+super_base + v_space*4 + base_height*4.75 +0.2,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p6 <- plotSignal(data = df_h3k27ac,
  params = c(chipseq_params, h3k27ac_range),
  fill = "#253494", linecolor = "#253494",
  x = 0.5, y =1.3+ super_base + v_space*4 + base_height*4.75 + 0.2 +0.2,
  height = base_height*0.75, width = figure_width)

annoYaxis(plot = p6, at = c(0, 40), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "H3K27ac", x = 0.1,
  y = 1.3+super_base + v_space*4 + base_height*4.75 + 0.2 + (base_height*0.75) / 2 +0.2,
  fontsize = annot_font_size, fontface = "bold",
  rot = 90,
  default.units = "inches")

plotText(label = "-log10(P)", x = 0.2,
  y = 1.3+super_base + v_space*4 + base_height*4.75 + 0.2 + (base_height*0.75) / 2 +0.2,
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
  x = 0.5, y = 1.3+super_base + v_space*5 + base_height*5.5 + 0.2 +0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "PBMC", x = 0.2,
  y = 1.3+super_base + v_space*5 + base_height*5.5 + 0.2 + (base_height*0.5)/2 +0.2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p8 <- plotRanges(
  data = t_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = 1.3+super_base + v_space*6 + base_height*6 + 0.2 +0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Naive T", x = 0.2,
  y = 1.3+super_base + v_space*6 + base_height*6 + 0.2 + (base_height*0.5)/2 +0.2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p9 <- plotRanges(
  data = b_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y =1.3+ super_base + v_space*7 + base_height*6.5 + 0.2 +0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Naive B", x = 0.2,
  y =1.3+ super_base + v_space*7 + base_height*6.5 + 0.2 + (base_height*0.5)/2 +0.2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p10 <- plotRanges(
  data = nk_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y =1.3+ super_base + v_space*8 + base_height*7 + 0.2 +0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "cNK", x = 0.2,
  y = 1.3+super_base + v_space*8 + base_height*7 + 0.2 + (base_height*0.5)/2 +0.2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p11 <- plotRanges(
  data = mono_atac, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes,
  assembly = "hg38",
  fill = colorby("type", palette = colorRampPalette(c("#253494", "lightgrey"))),
  linecolor = "fill",
  order = "random", collapse = TRUE,
  x = 0.5, y = 1.3+super_base + v_space*9 + base_height*7.5 + 0.2 +0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "Monocyte", x = 0.2,
  y = 1.3+super_base + v_space*9 + base_height*7.5 + 0.2 + (base_height*0.5)/2 +0.2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

annoGenomeLabel(plot= p8, x = 0.5,
  y = 1.3+super_base + v_space*9 + base_height*8.25 + 0.2 +0.2,
  fontsize = label_font_size, scale = "Mb", default.units = "inches")

pageGuideHide()

dev.off()

