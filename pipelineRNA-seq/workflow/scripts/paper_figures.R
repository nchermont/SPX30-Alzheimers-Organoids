#!/usr/bin/env Rscript
# Tables and figures of the paper (GROUND vs ISS analysis of the AD lines APP,
# PSEN1, PSEN2 and the wild-type line WTAD) from the STAR gene counts.
# Run by rule `paper_figures` (workflow/rules/paper.smk), or by hand from the
# repository root:
#   Rscript workflow/scripts/paper_figures.R [--out results/paper_figures]
#       [--counts results/counts/all.tsv] [--samples config/samples.tsv]
#       [--annotation data/gene_annotation.tsv]
#
# Colours: blue = up, orange = down; one colour per genotype and environment in
# every figure (docs/paper_figures.md).
#
# Software: workflow/envs/paper_figures.yaml

suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log)
  sink(log, type = "message")
  counts_file <- snakemake@input[["counts"]]
  samples_file <- snakemake@input[["samples"]]
  annotation_file <- snakemake@input[["annotation"]]
  out_dir <- snakemake@output[[1]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- function(name, default) {
    i <- match(paste0("--", name), args)
    if (is.na(i)) default else args[i + 1]
  }
  counts_file <- opt("counts", "results/counts/all.tsv")
  samples_file <- opt("samples", "config/samples.tsv")
  annotation_file <- opt("annotation", "data/gene_annotation.tsv")
  out_dir <- opt("out", "results/paper_figures")
}

PADJ <- 0.05  # significance threshold used throughout the paper
LFC <- 1      # |log2 fold change| threshold used throughout the paper
GENOTYPES <- c("APP", "PSEN1", "PSEN2")

out <- function(...) {
  path <- file.path(out_dir, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


##### colours ##################################################################

# blue = up (higher in ISS, or in AD / the AD line), orange = down; lighter
# tones for a second up or down set in Venn diagrams
geno <- c(WTAD = "#7F7F7F", APP = "#E69F00", PSEN1 = "#009E73", PSEN2 = "#CC79A7")
env <- c(GROUND = "#8C6D31", ISS = "#56B4E9")
up <- "blue"
down <- "orange"
up_light <- "#7F7FFF"
down_light <- "#FFD27F"
pal <- list(
  genotype = geno,
  heatmap_annotation = list(environment = env, condition = geno,
                            Highlight = c(Yes = "black", No = "white")),
  heatmap = colorRampPalette(c(down, "white", up))(100),
  volcano = c(Upregulated = up, Downregulated = down, `No change` = "grey"),
  volcano_all = c(`ISS up` = up, `GROUND up` = down, `No change` = "grey"),
  enrichment = c(low = down, high = up)
)

# genes marked in the highlighted interaction heatmaps and volcano plots
HIGHLIGHT <- list(
  APP = c("BAMBI", "SPATA6", "ACHE", "GRIA4", "LGI2", "RASGRF1", "GABRG2", "GAD2", "SLC32A1", "GAD1"),
  PSEN1 = c("SLC32A1", "GAD2", "GRIA4", "GABRG2", "GAD1", "SNAP25", "IGFBP2", "BAMBI"),
  PSEN2 = c("GAD1", "LBX1", "SLC32A1", "LHX5", "LHX1", "GAD2", "SNAP25")
)


##### data #####################################################################

counts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
sheet <- read.delim(samples_file, check.names = FALSE)
stopifnot(setequal(sheet$sample_name, colnames(counts)))
# sample-sheet order (the workflow's count matrix is sorted by sample name)
counts <- counts[, sheet$sample_name]

annotation <- read.delim(annotation_file, colClasses = "character",
                         na.strings = "NA", row.names = 1)
# biomaRt name, NA when biomaRt did not return the gene ("" = returned unnamed)
symbol_biomart <- function(ids) annotation[ids, "symbol_biomart"]
# name with the Ensembl id of the gene itself as fallback
named <- function(ids, symbols) ifelse(is.na(symbols) | symbols == "", ids, symbols)
# volcano labels: gene symbol (org.Hs.eg.db, else biomaRt); NA for unnamed genes,
# which are not labelled
label_symbol <- function(ids) {
  s <- annotation[ids, "symbol_orgdb"]
  b <- annotation[ids, "symbol_biomart"]
  s[is.na(s)] <- b[is.na(s)]
  s[!is.na(s) & s == ""] <- NA
  s
}

coldata <- function(samples) {
  s <- sheet[match(samples, sheet$sample_name), ]
  data.frame(
    row.names = samples,
    condition = factor(s$GENOTYPE, levels = c("WTAD", GENOTYPES)),
    environment = factor(s$ENVIRONMENT, levels = c("GROUND", "ISS")),
    disease = factor(ifelse(s$GENOTYPE == "WTAD", "WTAD", "AD"), levels = c("WTAD", "AD"))
  )
}

# GROUND and ISS samples, in sample-sheet order
gi <- colnames(counts)[sheet$ENVIRONMENT[match(colnames(counts), sheet$sample_name)] %in% c("GROUND", "ISS")]

dds_fit <- function(samples, design, min_count = 10, fit = TRUE) {
  x <- counts[, samples]
  if (!is.na(min_count)) x <- x[rowSums(x) > min_count, ]
  dds <- DESeqDataSetFromMatrix(x, coldata(samples), design)
  if (fit) DESeq(dds, quiet = TRUE) else dds
}

# condition-major (WTAD GROUND, WTAD ISS, APP GROUND, ...) and
# environment-major (GROUND WTAD, GROUND APP, ..., ISS WTAD, ...) column orders
cond_major <- function(cd) rownames(cd)[order(cd$condition, cd$environment)]
env_major <- function(cd) rownames(cd)[order(cd$environment, cd$condition)]

# sample names: full sheet name, the numbered name of vst_transformed_counts_all.csv
# and heatmap_app.pdf, and a standardised name keeping replicate and lane
cd_gi <- coldata(gi)
group <- paste(cd_gi$environment, cd_gi$condition, sep = "_")
sample_names <- data.frame(
  sample_name = gi,
  paper_short_name = paste(group, ave(seq_along(gi), group, FUN = seq_along), sep = "_"),
  standard_name = sub("^([A-Z]+_[A-Z0-9]+_[A-C])_S[0-9]+_(L00[0-9])$", "\\1_\\2", gi),
  environment = cd_gi$environment,
  genotype = cd_gi$condition,
  replicate = sub("^[A-Z]+_[A-Z0-9]+_([A-C])_.*", "\\1", gi),
  lane = sub(".*_(L00[0-9])$", "\\1", gi)
)
write.table(sample_names, out("sample_names.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
short_name <- setNames(sample_names$paper_short_name, sample_names$sample_name)


##### models ###################################################################

message("fitting DESeq2 models")
# genotype x environment (figures, interaction tables, ALL ISS vs GROUND)
dds_gi <- dds_fit(gi, ~ condition * environment)
# AD vs WTAD ignoring the environment (ALL AD vs WTAD)
dds_disease <- dds_fit(gi, ~ disease)
# ISS vs GROUND within each group, without gene pre-filtering
group_samples <- list(
  AD = gi[cd_gi$disease == "AD"], WTAD = gi[cd_gi$condition == "WTAD"],
  APP = gi[cd_gi$condition == "APP"], PSEN1 = gi[cd_gi$condition == "PSEN1"],
  PSEN2 = gi[cd_gi$condition == "PSEN2"]
)
dds_group <- lapply(group_samples, dds_fit, design = ~ environment, min_count = NA)

vst_blind <- assay(vst(dds_gi, blind = TRUE))    # PCA, top-1000 heatmap
vst_design <- assay(vst(dds_gi, blind = FALSE))  # interaction heatmaps


##### differential expression tables ###########################################

message("writing tables")
interaction_table <- function(dds, g) {
  r <- results(dds, name = sprintf("condition%s.environmentISS", g))
  r <- data.frame(ensembl_gene_id = rownames(r), as.data.frame(r))
  r <- r[order(r$ensembl_gene_id), ]
  r$external_gene_name <- symbol_biomart(r$ensembl_gene_id)
  rownames(r) <- NULL
  r
}
interaction <- list()
for (g in GENOTYPES) {
  interaction[[g]] <- interaction_table(dds_gi, g)
  write.csv(interaction[[g]], out(sprintf("interaction_%s_vs_WTAD_ground_vs_iss.csv", g)))
}

# per-group ISS vs GROUND; genes filtered by independent filtering get padj = 1
group_dir <- c(AD = "AD-WTAD", WTAD = "AD-WTAD", APP = "APP", PSEN1 = "PSEN1", PSEN2 = "PSEN2")
# the wild-type line is also written as the control group
CONTROL <- "WTAD"
group_res <- list()
summary_rows <- list()
for (g in names(dds_group)) {
  r <- as.data.frame(results(dds_group[[g]], contrast = c("environment", "ISS", "GROUND")))
  r <- r[r$baseMean > 0, ]
  summary_rows[[g]] <- data.frame(
    contrast = paste0(g, "-ISS_vs_GROUND"),
    n_tested = sum(!is.na(r$padj)),
    n_deg_padj = sum(r$padj < PADJ, na.rm = TRUE),
    n_deg_padj_lfc = sum(r$padj < PADJ & abs(r$log2FoldChange) > LFC, na.rm = TRUE),
    n_up = sum(r$padj < PADJ & r$log2FoldChange > LFC, na.rm = TRUE),
    n_down = sum(r$padj < PADJ & r$log2FoldChange < -LFC, na.rm = TRUE)
  )
  r$padj[is.na(r$padj)] <- 1
  r$gene <- rownames(r)
  r$ensembl <- rownames(r)
  r$symbol <- named(r$ensembl, annotation[r$ensembl, "symbol_orgdb"])
  group_res[[g]] <- r
  write.csv(r, out(group_dir[[g]], sprintf("%s_iss_vs_ground_all.csv", tolower(g))), row.names = FALSE)
  if (g == CONTROL) write.csv(r, out("control", "control_iss_vs_ground_all.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, summary_rows), out("DEG_summary_individual_compare.csv"), row.names = FALSE)

ordered_table <- function(r, up, down) {
  r <- as.data.frame(r)
  r <- r[!is.na(r$padj), ]
  r <- r[order(r$padj), ]
  data.frame(
    ensembl = rownames(r),
    gene_symbol = named(rownames(r), annotation[rownames(r), "symbol_all"]),
    r[, c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
    status = ifelse(r$padj < PADJ & r$log2FoldChange > LFC, up,
                    ifelse(r$padj < PADJ & r$log2FoldChange < -LFC, down, "No change"))
  )
}
# NB: in ~condition * environment this coefficient is the ISS effect in WTAD
all_iss <- ordered_table(results(dds_gi, name = "environment_ISS_vs_GROUND"), "ISS up", "GROUND up")
all_ad <- ordered_table(results(dds_disease, contrast = c("disease", "AD", "WTAD")), "AD up", "WTAD up")
write.csv(all_iss, out("ALL", "DESeq2_ISS_vs_GROUND_results_ordered.csv"), row.names = FALSE)
write.csv(all_ad, out("ALL", "DESeq2_AD_vs_WTAD_results_ordered.csv"), row.names = FALSE)

write_genes <- function(genes, ...) writeLines(genes, out(...))
up_in <- function(t, s) t$gene_symbol[t$status == s]
write_genes(up_in(all_iss, "ISS up"), "ALL", "genes_up_in_ISS.txt")
write_genes(up_in(all_iss, "GROUND up"), "ALL", "genes_up_in_GROUND.txt")
write_genes(up_in(all_ad, "AD up"), "ALL", "genes_up_in_AD.txt")
write_genes(up_in(all_ad, "WTAD up"), "ALL", "genes_up_in_WTAD.txt")

deg_symbols <- function(r, sign) unique(r$symbol[r$padj < PADJ & sign * r$log2FoldChange > LFC])
iss_lists <- list(
  AD_Up = deg_symbols(group_res$AD, 1), WTAD_Up = deg_symbols(group_res$WTAD, 1),
  WTAD_Down = deg_symbols(group_res$WTAD, -1), AD_Down = deg_symbols(group_res$AD, -1)
)
write_genes(iss_lists$AD_Up, "AD-WTAD", "AD_genes_up_in_ISS.txt")
write_genes(iss_lists$AD_Down, "AD-WTAD", "AD_genes_down_in_ISS.txt")
write_genes(iss_lists$WTAD_Up, "AD-WTAD", "WTAD_genes_up_in_ISS.txt")
write_genes(iss_lists$WTAD_Down, "AD-WTAD", "WTAD_genes_down_in_ISS.txt")
write.table(data.frame(Gene = intersect(iss_lists$AD_Up, iss_lists$WTAD_Up)),
            out("AD-WTAD", "AD_WTAD_up_intersection.txt"), row.names = FALSE)
write.table(data.frame(Gene = intersect(iss_lists$AD_Down, iss_lists$WTAD_Down)),
            out("AD-WTAD", "AD_WTAD_down_intersection.txt"), row.names = FALSE)


##### normalised and transformed counts ########################################

norm <- counts(dds_gi, normalized = TRUE)
keep <- sort(rownames(norm)[!is.na(symbol_biomart(rownames(norm)))])
write.csv(data.frame(ensembl_gene_id = keep, external_gene_name = symbol_biomart(keep),
                     norm[keep, ], check.names = FALSE),
          out("normalized_counts_all_genes.csv"), row.names = FALSE)

zscore <- function(m) t(scale(t(m)))
top1000 <- rownames(vst_blind)[order(rowVars(vst_blind), decreasing = TRUE)[1:1000]]
z_top <- zscore(vst_blind[top1000, ])
write.csv(z_top[, env_major(cd_gi)], out("heatmap_top1000_variable_genes.csv"))
# row z-scores of the VST values (file name kept from the first version of the paper)
top_symbol <- ifelse(is.na(symbol_biomart(top1000)), top1000, symbol_biomart(top1000))
write.csv(data.frame(GeneSymbol = top_symbol, z_top[, cond_major(cd_gi)], check.names = FALSE),
          out("vst_transformed_counts_top1000_genes.csv"))
write.csv(data.frame(GeneSymbol = top_symbol, norm[top1000, cond_major(cd_gi)], check.names = FALSE),
          out("normalized_counts_top1000_genes.csv"))
z_all <- zscore(vst_blind)[, cond_major(cd_gi)]
colnames(z_all) <- short_name[colnames(z_all)]
write.csv(data.frame(`Sample ID` = rownames(z_all), z_all, check.names = FALSE),
          out("vst_transformed_counts_all.csv"), row.names = FALSE, quote = FALSE)


##### PCA ######################################################################

message("plotting")
pca_plot <- function(v, cd, file, labels = FALSE, legend_shape_first = TRUE) {
  top <- order(rowVars(v), decreasing = TRUE)[1:min(500, nrow(v))]
  pca <- prcomp(t(v[top, ]))
  var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cd[colnames(v), ],
                  label = paste(cd[colnames(v), "condition"], cd[colnames(v), "environment"], sep = ":"))
  d$condition <- droplevels(d$condition)
  p <- ggplot(d, aes(PC1, PC2, colour = condition, shape = environment)) +
    geom_point(size = if (legend_shape_first) 4 else 3, alpha = 0.8) +
    labs(x = sprintf("PC1: %d%% variance", var[1]), y = sprintf("PC2: %d%% variance", var[2]),
         title = "PCA: Condition + Environment", colour = NULL, shape = NULL) +
    theme_bw(base_size = if (legend_shape_first) 18 else 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"))
  p <- p + guides(shape = guide_legend(order = if (legend_shape_first) 1 else 2),
                  colour = guide_legend(order = if (legend_shape_first) 2 else 1))
  p <- p + scale_colour_manual(values = pal$genotype)
  if (labels) p <- p + geom_text_repel(aes(label = label), size = 3.5, show.legend = FALSE, max.overlaps = Inf)
  ggsave(out(file), p, width = 10, height = 6, dpi = 300, bg = "white")
}
pca_plot(vst_blind, cd_gi, "pca_plot.png")
pca_plot(vst_blind, cd_gi, "pca_plot_with_labels.png", labels = TRUE)
for (g in GENOTYPES) {
  pair <- gi[cd_gi$condition %in% c("WTAD", g)]
  v <- assay(vst(dds_fit(pair, ~ condition * environment, fit = FALSE), blind = TRUE))
  pca_plot(v, cd_gi, sprintf("pca_plot_%s.png", g), legend_shape_first = FALSE)
}


##### heatmaps #################################################################

# pheatmap stacks annotation tracks bottom-up: environment is drawn on top
heatmap_annotation <- function(samples) {
  a <- cd_gi[samples, c("condition", "environment")]
  a$condition <- droplevels(a$condition)
  a
}
annotation_colours <- function(a, extra = NULL) {
  c(lapply(setNames(names(a), names(a)), function(n) pal$heatmap_annotation[[n]][levels(a[[n]])]), extra)
}

png(out("heatmap_top1000_variable_genes.png"), width = 7, height = 7, units = "in", res = 300)
a <- heatmap_annotation(env_major(cd_gi))
pheatmap(z_top[, env_major(cd_gi)], color = pal$heatmap, cluster_cols = FALSE, show_rownames = FALSE,
         annotation_col = a, annotation_colors = annotation_colours(a),
         main = "Top 1000 Most Variable Genes", fontsize = 12, fontsize_col = 8, border_color = NA)
dev.off()

centred <- vst_design - rowMeans(vst_design)
for (g in GENOTYPES) {
  t <- interaction[[g]]
  t <- t[order(t$padj), ]
  t <- t[!is.na(t$padj) & !is.na(t$external_gene_name) & t$external_gene_name != "", ]
  top <- head(t, 30)
  samples <- env_major(cd_gi[cd_gi$condition %in% c("WTAD", g), ])
  values <- centred[top$ensembl_gene_id, samples]
  rownames(values) <- top$external_gene_name
  write.csv(values, out(sprintf("heatmap_interaction_%s_ground_vs_iss_values.csv", g)))

  a <- heatmap_annotation(samples)
  title <- sprintf("Top Interaction Genes: %s vs WTAD × ISS", g)
  png(out(sprintf("heatmap_interaction_%s_ground_vs_iss.png", g)), width = 7, height = 7, units = "in", res = 300)
  pheatmap(values, color = pal$heatmap, cluster_cols = FALSE, annotation_col = a,
           annotation_colors = annotation_colours(a), main = title, fontsize = 12, fontsize_row = 8)
  dev.off()

  hl <- rownames(values) %in% HIGHLIGHT[[g]]
  row_a <- data.frame(Highlight = factor(ifelse(hl, "Yes", "No"), levels = c("Yes", "No")), row.names = rownames(values))
  png(out(sprintf("heatmap_interaction_%s_ground_vs_iss_highlighted.png", g)), width = 7, height = 7, units = "in", res = 300)
  pheatmap(values, color = pal$heatmap, cluster_cols = FALSE, annotation_col = a, annotation_row = row_a,
           annotation_colors = annotation_colours(a, list(Highlight = pal$heatmap_annotation$Highlight)),
           labels_row = ifelse(hl, paste0("*", rownames(values), "*"), rownames(values)), main = title,
           fontsize = 12, fontsize_row = 8)
  dev.off()
}

# APP heatmap with numbered sample names, columns clustered (heatmap_app.pdf)
values <- read.csv(out("heatmap_interaction_APP_ground_vs_iss_values.csv"), row.names = 1, check.names = FALSE)
group_a <- data.frame(group = factor(paste(cd_gi[colnames(values), "environment"], cd_gi[colnames(values), "condition"], sep = "_")),
                      row.names = short_name[colnames(values)])
colnames(values) <- short_name[colnames(values)]
group_colours <- c(GROUND_APP = "#F0C566", GROUND_WTAD = "#BFBFBF", ISS_APP = "#B37000", ISS_WTAD = "#4D4D4D")
pdf(out("heatmap_app.pdf"), width = 7, height = 6)
pheatmap(values, color = pal$heatmap, annotation_col = group_a, annotation_colors = list(group = group_colours),
         annotation_names_col = FALSE, fontsize_row = 7, fontsize_col = 7, border_color = NA)
dev.off()


##### volcano plots ############################################################

regulation <- function(r) {
  factor(ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FoldChange > LFC, "Upregulated",
                ifelse(!is.na(r$padj) & r$padj < PADJ & r$log2FoldChange < -LFC, "Downregulated", "No change")),
         levels = c("Downregulated", "No change", "Upregulated"))
}

volcano_interaction <- function(r, g, file, highlight = NULL, cap = FALSE) {
  d <- data.frame(x = r$log2FoldChange, y = -log10(r$padj), symbol = r$external_gene_name, threshold = regulation(r))
  d <- d[!is.na(d$y), ]
  xlab <- "log2FoldChange"
  ylab <- "-log10(padj)"
  if (cap) {
    d$x <- pmax(pmin(d$x, 5), -5)
    d$y <- pmin(d$y, 50)
    xlab <- "log2FoldChange_capped"
    ylab <- "neglog10padj_capped"
  }
  p <- ggplot(d, aes(x, y, colour = threshold)) +
    geom_point(alpha = 0.5, size = 2) +
    scale_colour_manual(values = pal$volcano) +
    labs(title = sprintf("%s vs WTAD interaction with ISS", g), x = xlab, y = ylab) +
    theme_minimal(base_size = 14)
  if (!is.null(highlight)) {
    h <- d[d$symbol %in% highlight, ]
    p <- p + geom_text_repel(data = h, aes(label = symbol), size = 6, show.legend = FALSE, max.overlaps = Inf)
  }
  ggsave(out(file), p, width = 6.67, height = 6.67, dpi = 300, bg = "white")
}
for (g in GENOTYPES) {
  r <- interaction[[g]]
  volcano_interaction(r, g, sprintf("volcano_interaction_%s_ground_vs_iss.png", g))
  volcano_interaction(r, g, sprintf("volcano_interaction_%s_ground_vs_iss_highlighted.png", g), HIGHLIGHT[[g]])
  volcano_interaction(r, g, sprintf("volcano_interaction_%s_ground_vs_iss_highlighted_capped.png", g), HIGHLIGHT[[g]], cap = TRUE)
}

# per-group ISS vs GROUND: log2FC capped at +-5, -log10 FDR at 50, top 10 up/down labelled
volcano_group <- function(r, title, file) {
  d <- data.frame(x = pmax(pmin(r$log2FoldChange, 5), -5), y = pmin(-log10(r$padj), 50),
                  padj = r$padj, symbol = label_symbol(r$ensembl), threshold = regulation(r))
  # top 10 by FDR among named genes, each symbol labelled once
  top10 <- function(k) {
    x <- d[k & !is.na(d$symbol), ]
    x <- x[order(x$padj), ]
    head(x[!duplicated(x$symbol), ], 10)
  }
  lab <- rbind(top10(d$threshold == "Upregulated"), top10(d$threshold == "Downregulated"))
  p <- ggplot(d, aes(x, y, colour = threshold)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_text_repel(data = lab, aes(label = symbol), size = 3.5, show.legend = FALSE, max.overlaps = Inf) +
    scale_colour_manual(values = pal$volcano) +
    expand_limits(y = max(d$y, na.rm = TRUE) + 10) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted p-value (FDR)") +
    theme_minimal(base_size = 13)
  ggsave(out(file), p, width = 8, height = 6, dpi = 300, bg = "white")
}
for (g in names(group_res)) {
  volcano_group(group_res[[g]], sprintf("%s: ISS vs GROUND", g),
                file.path(group_dir[[g]], sprintf("volcano_%s_iss_vs_ground.png", tolower(g))))
}
volcano_group(group_res[[CONTROL]], sprintf("Control (%s): ISS vs GROUND", CONTROL),
              file.path("control", "volcano_control_iss_vs_ground.png"))

d <- all_iss
d$y <- -log10(d$padj)
d$label <- ifelse(d$gene_symbol == d$ensembl, label_symbol(d$ensembl), d$gene_symbol)
top10 <- function(x) head(x[!is.na(x$label) & !duplicated(x$label), ], 10)  # d is ordered by padj
lab <- rbind(top10(d[d$status == "ISS up", ]), top10(d[d$status == "GROUND up", ]))
p <- ggplot(d, aes(log2FoldChange, y, colour = status)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed") +
  geom_text_repel(data = lab, aes(label = label), size = 4, show.legend = FALSE, max.overlaps = Inf) +
  scale_colour_manual(values = pal$volcano_all) +
  labs(title = "Volcano plot: ISS vs GROUND", x = "log2 Fold Change (ISS / GROUND)",
       y = "-log10 adjusted p-value", colour = "Regulation") +
  theme_minimal(base_size = 14)
ggsave(out("ALL", "volcano_ISS_vs_GROUND.png"), p, width = 16, height = 9, dpi = 300, bg = "white")


##### Venn diagrams ############################################################

venn <- function(sets, colours, file, main = NULL, sub = NULL, ...) {
  venn.diagram(sets, filename = out(file), imagetype = "png", height = 8, width = 8,
               units = "in", resolution = 300, fill = unname(colours), alpha = 0.6,
               main = main, sub = sub, disable.logging = TRUE, ...)
}
iss_up <- up_in(all_iss, "ISS up")
ground_up <- up_in(all_iss, "GROUND up")
ad_up <- up_in(all_ad, "AD up")
wtad_up <- up_in(all_ad, "WTAD up")
venn(list(ISS = iss_up, GROUND = ground_up), c(up, down), file.path("ALL", "venn_ISS_vs_GROUND.png"),
     "Differentially Expressed Genes", "ISS vs GROUND")
venn(list(AD = ad_up, WTAD = wtad_up), c(up, down), file.path("ALL", "venn_AD_vs_WTAD.png"),
     "Differentially Expressed Genes", "AD vs WTAD")
three <- list(cat.pos = c(-40, 0, 40), cat.dist = c(0.05, 0.05, 0.05))
do.call(venn, c(list(list(GROUND_up = ground_up, AD_up = ad_up, ISS_up = iss_up), c(down, up_light, up),
                      file.path("ALL", "venn_ISS_GROUND_ADUP.png"), "Differentially Expressed Genes",
                      "ISS/GROUND and AD/WTAD"), three))
do.call(venn, c(list(list(GROUND_up = ground_up, AD_down = wtad_up, ISS_up = iss_up), c(down, down_light, up),
                      file.path("ALL", "venn_ISS_GROUND_ADDOWN.png"), "Differentially Expressed Genes",
                      "ISS/GROUND and AD/WTAD"), three))
venn(list(AD = iss_lists$AD_Up, WTAD = iss_lists$WTAD_Up), c(up, up_light),
     file.path("AD-WTAD", "venn_up_ad_vs_wtad.png"), "Upregulated: AD vs WTAD", cex = 2, cat.cex = 2, main.cex = 2)
venn(list(AD = iss_lists$AD_Down, WTAD = iss_lists$WTAD_Down), c(down, down_light),
     file.path("AD-WTAD", "venn_down_ad_vs_wtad.png"), "Downregulated: AD vs WTAD", cex = 2, cat.cex = 2, main.cex = 2)
# VennDiagram places 4 sets as 1, 3, 4, 2 from left to right
venn(iss_lists[c("AD_Up", "AD_Down", "WTAD_Up", "WTAD_Down")],
     c(up, down, up_light, down_light),
     file.path("AD-WTAD", "venn_ad_wtad_up_down.png"), "AD/WTAD Up/Down Overlap", cex = 1.5, cat.cex = 1.5, main.cex = 1.5)


##### GO enrichment ############################################################

message("GO enrichment")
for (g in GENOTYPES) {
  t <- interaction[[g]]
  universe <- unique(t$external_gene_name[!is.na(t$external_gene_name) & t$external_gene_name != ""])
  sig <- t$external_gene_name[!is.na(t$padj) & t$padj < PADJ & abs(t$log2FoldChange) > LFC]
  sig <- unique(sig[!is.na(sig) & sig != ""])
  ego <- enrichGO(sig, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", universe = universe,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  write.csv(as.data.frame(ego), out(sprintf("%s_vs_WTAD_ISS_enrichGO_BP.csv", g)), row.names = FALSE)
  fill <- scale_fill_gradient(low = pal$enrichment[["low"]], high = pal$enrichment[["high"]])
  colour <- scale_colour_gradient(low = pal$enrichment[["low"]], high = pal$enrichment[["high"]])
  ggsave(out(sprintf("%s_vs_WTAD_ISS_dotplot_BP.png", g)),
         suppressMessages(dotplot(ego, showCategory = 20) + fill + colour + ggtitle("GO BP Dotplot")),
         width = 10, height = 6.67, dpi = 300, bg = "white")
  ggsave(out(sprintf("%s_vs_WTAD_ISS_barplot_BP.png", g)),
         suppressMessages(barplot(ego, showCategory = 20) + fill + ggtitle("GO BP Barplot")),
         width = 12, height = 8, dpi = 300, bg = "white")
}

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("done: ", out_dir)
