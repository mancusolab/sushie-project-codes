# this script is to create main figure 3
library(tidyverse)
library(glue)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(RColorBrewer)

rnaseq_meta <- read_tsv("~/Documents/github/data/sushie_results/real2/rnaseq_normal.meta_cs.tsv.gz")

# case study URGCP
target_gene <- "ENSG00000203875_SNHG5"
target_snp <- "rs1059307"

# check if this gene has credible sets in other results

# zgrep "SNHG5" ~/Documents/github/data/sushie_results/real2/rnaseq_susiex_cs.tsv.gz
# zgrep "SNHG5" ~/Documents/github/data/sushie_results/real2/rnaseq_mesusie_cs.tsv.gz
# zgrep "SNHG5" ~/Documents/github/data/sushie_results/real2/rnaseq_xmap_cs.tsv.gz

df_eqtl_eur <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.EUR.qtlscan.PHENO1.glm.linear"))

df_eqtl_afr <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.AFR.qtlscan.PHENO1.glm.linear"))

df_eqtl_his <- read_tsv(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.HIS.qtlscan.PHENO1.glm.linear"))

dd <- df_eqtl_eur %>%
  dplyr::select(CHROM = `#CHROM`, POS, ID, REF, ALT,
    EUR_BETA = BETA, EUR_SE = SE, EUR_T_STAT = T_STAT) %>%
  left_join(df_eqtl_afr %>%
      dplyr::select(ID, AFR_BETA = BETA, AFR_SE = SE, AFR_T_STAT = T_STAT)) %>%
  left_join(df_eqtl_his %>%
      dplyr::select(ID, HIS_BETA = BETA, HIS_SE = SE, HIS_T_STAT = T_STAT))

# write_tsv(dd, "./tables/s11.tsv")

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

target_chrom <- "chr6"
# this is flank window, not actual tss and tes
target_tss <- 85150491
target_tes <- 86178931
snp_pos <- 85678170

df_h3k27ac <- readBigwig(glue(
  "~/Documents/github/data/sushie_results/case/{target_gene}/SNHG5_h3k27ac.bigWig"),
  chrom = target_chrom, chromstart = target_tss, chromend = target_tes)

atac_seq <- read_tsv(
  glue("~/Documents/github/data/sushie_results/case/snATAC-seq-frozen-peaks_grch38.bed"), col_names = FALSE) %>%
  dplyr::select(chrom = X1, start = X2, end = X3) %>%
  mutate(type = factor(ifelse(chrom == target_chrom & start <= snp_pos &
      snp_pos <= end, "Yes", "No"), levels = c("Yes", "No")))

chipseq_params <- pgParams(chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, 
  assembly = "hg38", x = 3.5, width = 1.5, default.units = "inches")

h3k27ac_range <- pgParams(range = c(min(df_h3k27ac$score),
  max(df_h3k27ac$score)), assembly = "hg38")

df_ld1 <- read_tsv(
  glue("~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.0.trueLD.tsv")) %>%
  filter(counts != 1)

df_ld2 <- read_tsv(
  glue("~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.1.trueLD.tsv")) %>%
  filter(counts != 1) 

df_ld3 <- read_tsv(
  glue("~/Documents/github/data/sushie_results/case/{target_gene}/{target_gene}.2.trueLD.tsv")) %>%
  filter(counts != 1) 

max_counts <- max(df_ld1$counts, df_ld2$counts, df_ld3$counts)

df_pairs <- tibble("chrom1" = target_chrom, "start1" = snp_pos, "end1" = snp_pos,
  "chrom2" = target_chrom, "start2" = snp_pos, "end2" = snp_pos)

png(filename="./plots/s16.png", width = 7, height = 5.4, units="in", res=300)

base_height <- 0.5
v_space <- 0.1
super_base <- 0.2
label_font_size <- 8
axis_font_size <- 7
legend_font_size <- 11
figure_width <- 6.2
annot_font_size <- 6

pageCreate(width = 7, height = 5.4, default.units = "inches")

plotText(label = "A", x = 0.1, y = 0.1,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p_ld1 <- plotHicTriangle(df_ld1, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0, max_counts),
  x=0.5, y=super_base + 0.3, width = 1.8,
  height = base_height*1.8, just = c("left"),
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld1, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "EUR", x = 0.5+1.8/2, y =super_base + 0.85,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

p_ld2<- plotHicTriangle(df_ld2, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0,max_counts),
  x=0.5 + 1.8 + 0.2, y=super_base + 0.3, width = 1.8,
  height = base_height*1.8, just = "left",
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld2, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "AFR", x = 0.5 + 1.8 + 0.2 +1.8/2,
  y =super_base + 0.85,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

p_ld3 <- plotHicTriangle(df_ld3, chrom=target_chrom, resolution = 25000,
  chromstart = target_tss, chromend = target_tes, zrange = c(0,max_counts),
  x=0.5 + 1.8*2 + 0.2*2, y=super_base + 0.3, width = 1.8,
  height = base_height*1.8, just = "left",
  palette = colorRampPalette(brewer.pal(n = 9, "BuPu")),
  default.units = "inches")

pixels <- annoPixels(
  plot = p_ld3, data = df_pairs, type = "box",
  just="center", shift=0.5, col="red",
  lwd=5, lty=1
)

plotText(label = "HIS", x = 0.5 + 1.8*2 + 0.2*2 + 1.8/2,
  y =super_base + 0.85,
  fontsize = legend_font_size, fontface = "bold", just = "top",
  default.units = "inches")

annoHeatmapLegend(
  plot = p_ld3,
  x = 0.5 + 1.8*3 + 0.2*2 + 0.1, y = super_base + 0.3,
  width = 0.12, height = base_height*1.8,
  just = c("left"), default.units = "inches"
)

p1 <- plotManhattan(
  data = df_eur, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill = "grey", assembly = "hg38", x = 0.5, y = 1.1 +super_base, width = figure_width,
  height = base_height)

annoYaxis(plot = p1, axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "EUR", x = 6.8, y = 1.1 +super_base + base_height / 2,
  rot = 270, fontsize = label_font_size, fontface = "bold", default.units = "inches")

p2 <- plotManhattan(
  data = df_afr, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"), fill="grey",
  assembly = "hg38", x = 0.5, y = 1.1 +super_base + v_space + base_height,
  width = figure_width, height = base_height)

annoYaxis(plot = p2, axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "AFR", x = 6.8,
  y =1.1 + super_base + v_space + base_height + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

p3 <- plotManhattan(
  data = df_his, chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes,
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"), fill="grey",
  assembly = "hg38", x = 0.5,
  y = 1.1 +super_base + v_space*2 + base_height*2,
  width = figure_width, height = base_height
)

annoYaxis(plot = p3, axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "HIS", x = 6.8,
  y = 1.1 +super_base + v_space*2 + base_height*2 + base_height/2, rot = 270,
  fontsize = label_font_size, fontface = "bold", default.units = "inches")

plotText(label = "-log10(P)", x = 0.2,
  y = 1.1 +super_base + (v_space*2 + base_height*3) / 2, rot = 90,
  fontsize = annot_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p4 <- plotManhattan(data = df_weight, chrom = target_chrom,
  chromstart = target_tss, chromend = target_tes, trans = "",
  leadSNP = list(snp = "", pch = 18, cex = 0.6, fill = "red"),
  fill="grey", assembly = "hg38", x = 0.5,
  y = 1.1 +super_base + v_space*3 + base_height*3,
  width = figure_width, height = base_height*0.75)

annoYaxis(plot = p4, at = c(0, 1), axisLine = TRUE, fontsize = axis_font_size)

plotText(label = "PIP", x = 0.2,
  y = 1.1 +super_base + v_space*3 + base_height*3 + (base_height*0.75) / 2,
  rot = 90, fontsize = annot_font_size, fontface = "bold", default.units = "inches")

p5 <- plotGenes(assembly = "hg38", chrom = target_chrom, chromstart = target_tss,
  chromend = target_tes, geneHighlights = data.frame("gene" = c("SNHG5"),
    "color" = c("orange")),
  geneBackground = "grey", fontsize = axis_font_size, x = 0.5,
  y = 1.1 +super_base + v_space*4 + base_height*3.75,
  width=figure_width, height = base_height)

plotText(label = "B", x = 0.1, y = 1.1 +super_base + v_space*4 + base_height*4.75,
  fontsize = legend_font_size, fontface = "bold", just = "center",
  default.units = "inches")

p6 <- plotSignal(data = df_h3k27ac,
  params = c(chipseq_params, h3k27ac_range),
  fill = "#253494", linecolor = "#253494",
  x = 0.5, y = 1.1 +super_base + v_space*4 + base_height*4.75 + 0.2,
  height = base_height*0.75, width = figure_width)

annoYaxis(plot = p6, at = c(0, 25), axisLine = TRUE,
  fontsize = axis_font_size)

plotText(label = "H3K27ac", x = 0.1,
  y = 1.1 +super_base + v_space*4 + base_height*4.75 + 0.2 + (base_height*0.75) / 2,
  fontsize = annot_font_size, fontface = "bold",
  rot = 90,
  default.units = "inches")

plotText(label = "-log10(P)", x = 0.2,
  y = 1.1 +super_base + v_space*4 + base_height*4.75 + 0.2 + (base_height*0.75) / 2,
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
  x = 0.5, y = 1.1 +super_base + v_space*5 + base_height*5.5 + 0.2,
  width = figure_width, height = base_height*0.5,
  default.units = "inches"
)

plotText(label = "PBMC", x = 0.2,
  y = 1.1 +super_base + v_space*5 + base_height*5.5 + 0.2 + (base_height*0.5)/2,
  rot = 90,
  fontsize = annot_font_size, fontface = "bold", default.units = "inches")

annoGenomeLabel(plot= p7, x = 0.5,
  y = 1.1 +super_base + v_space*5 + base_height*6 + 0.3,
  fontsize = label_font_size, scale = "Mb", default.units = "inches")

pageGuideHide()

dev.off()

