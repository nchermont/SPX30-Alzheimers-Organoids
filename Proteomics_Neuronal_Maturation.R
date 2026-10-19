# ============================================================
# PROTEOMICS — NEURONAL DIFFERENTIATION / MATURATION ANALYSIS
# ============================================================

library(tidyverse)
library(openxlsx)
library(pheatmap)
library(ggrepel)
library(patchwork)
library(limma)

# =========================
# 1. SETTINGS
# =========================

proteomics_file <- "/mnt/data/Notch inhibition_log10_2(1).csv"
output_dir <- "Proteomics_Neuronal_Maturation2"

id_col <- "symbol"
already_log2 <- TRUE
max_missing_fraction <- 0.40

group_order <- c(
  "NDC2_DMSO", "NDC2_DAPT",
  "sAD1_DMSO", "sAD1_DAPT",
  "sAD2_DMSO", "sAD2_DAPT"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Tables"), showWarnings = FALSE)

# =========================
# 2. READ DATA
# =========================

prot_raw <- read.csv(proteomics_file, check.names = FALSE)

sample_cols <- grep("DMSO|DAPT", colnames(prot_raw), value = TRUE)

meta <- data.frame(
  sample = sample_cols,
  genotype = sub("_(DMSO|DAPT)_.*", "", sample_cols),
  treatment = sub(".*_(DMSO|DAPT)_.*", "\\1", sample_cols)
)

meta$group <- paste(meta$genotype, meta$treatment, sep = "_")
meta$group <- factor(meta$group, levels = group_order)
meta$genotype <- factor(meta$genotype, levels = c("NDC2", "sAD1", "sAD2"))
meta$treatment <- factor(meta$treatment, levels = c("DMSO", "DAPT"))

# Colors used across plots
group_colors <- c(
  "DMSO" = "#80B1D3",
  "DAPT" = "#C13C91"
)


# =========================
# 3. PREPARE MATRIX
# =========================

id_col_clean <- id_col

prot <- prot_raw %>%
  dplyr::filter(!is.na(.data[[id_col_clean]])) %>%
  dplyr::mutate(Gene = as.character(.data[[id_col_clean]])) %>%
  dplyr::select(Gene, dplyr::all_of(sample_cols))

prot_mat <- prot %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_cols), ~ median(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  )

mat <- prot_mat %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

mode(mat) <- "numeric"

if (!already_log2) {
  mat <- log2(mat + 1)
}

keep <- rowMeans(is.na(mat)) <= max_missing_fraction
mat <- mat[keep, ]

impute_min <- function(x) {
  if (all(is.na(x))) return(x)
  x[is.na(x)] <- min(x, na.rm = TRUE) - 1
  return(x)
}

mat_imp <- t(apply(mat, 1, impute_min))

mat_z <- t(scale(t(mat_imp)))
mat_z[is.na(mat_z)] <- 0

mat_imp <- mat_imp[, meta$sample, drop = FALSE]
mat_z <- mat_z[, meta$sample, drop = FALSE]

# =========================
# 4. CELL-TYPE / MATURATION MODULES
# =========================

modules <- list(
  Radial_glia = c(
    "NES","HES1","HES5","NOTCH1","PROM1","SOX2","PAX6",
    "VIM","FABP7","HOPX","PTPRZ1","TNC","FAM107A","LIFR"
  ),
  
  Cycling_progenitors = c(
    "MKI67","TOP2A","MCM2","PCNA"
  ),
  
  Intermediate_progenitors = c(
    "EOMES","NEUROD1","NEUROD2","NEUROG1","NEUROG2",
    "ASCL1","BTG2","DLL3","NHLH1","NHLH2","PPP1R17"
  ),
  
  Excitatory_neurons = c(
    "SLC17A6","SLC17A7","TBR1","SATB2","BCL11B",
    "CAMK2A","NPTX1","NRN1","NRXN1","GRIN1",
    "GRIN2B","GRIA1","GRIA2","FEZF2","SOX5","FOXP2",
    "CUX1","CUX2","POU3F2","POU3F3","RORB","LHX2"
  ),
  
  Inhibitory_interneurons = c(
    "GAD1","GAD2","SLC32A1","DLX1","DLX2","DLX5","DLX6",
    "ARX","SP8","PBX3","SST","CALB2","RELN","VIP",
    "LAMP5","LHX6","SOX6","PVALB","ERBB4","NR2F1","NR2F2"
  ),
  
  Astrocyte = c(
    "GFAP","AQP4","S100B","ALDH1L1","SLC1A2","SLC1A3",
    "SOX9","GJA1","GLUL","APOE","CLU","FGFR3","NDRG2",
    "ALDOC","SPARCL1"
  ),
  
  Reactive_astrocyte = c(
    "GFAP","VIM","CD44","C3","SERPINA3","CHI3L1",
    "STAT3","IFITM3","APOE","CLU","LCN2"
  ),
  
  Oligodendrocyte_lineage = c(
    "PDGFRA","CSPG4","OLIG1","OLIG2","SOX10","NKX2-2",
    "GPR17","VCAN","PTPRZ1","BCAN","MYRF","BCAS1",
    "ENPP6","TSPAN2","CNP","MBP","PLP1","MOG","MAG",
    "MOBP","CLDN11","MAL","ERMN","OPALIN"
  ),
  
  Choroid_plexus = c(
    "TTR","AQP1","KRT8","KRT18","KRT19","OTX2",
    "TPBG","RSPO2","RSPO3","CLDN2","IGFBP7"
  ),
  
  Mesenchymal_like = c(
    "COL1A1","COL1A2","COL3A1","COL6A1","COL6A2",
    "DCN","LUM","SPARC","POSTN","FN1","TAGLN",
    "ACTA2","ANXA2","SERPINH1"
  )
)
# =========================
# 4B. CALCULATE MODULE SCORES
# =========================

score_module <- function(mat_z, genes) {
  genes_present <- intersect(genes, rownames(mat_z))
  
  if (length(genes_present) < 2) {
    return(rep(NA_real_, ncol(mat_z)))
  }
  
  colMeans(mat_z[genes_present, , drop = FALSE], na.rm = TRUE)
}

module_scores <- purrr::map_dfc(names(modules), function(m) {
  tibble::tibble(!!m := score_module(mat_z, modules[[m]]))
}) %>%
  dplyr::mutate(sample = colnames(mat_z)) %>%
  dplyr::left_join(meta, by = "sample") %>%
  dplyr::relocate(sample, genotype, treatment, group)

write.xlsx(
  module_scores,
  file.path(output_dir, "Tables", "celltype_maturation_module_scores.xlsx")
)

module_detection <- purrr::imap_dfr(modules, function(genes, module) {
  tibble::tibble(
    module = module,
    gene = genes,
    detected = gene %in% rownames(mat_z)
  )
})

write.xlsx(
  module_detection,
  file.path(output_dir, "Tables", "celltype_module_gene_detection.xlsx")
)
# ============================================================
# 4C. SYNAPSE / EXCITATORY-INHIBITORY BALANCE MODULE
# ============================================================

synapse_modules <- list(
  
  Excitatory_identity = c(
    "SLC17A6","SLC17A7",
    "TBR1","SATB2","BCL11B",
    "CUX1","CUX2","RORB",
    "FEZF2","SOX5","FOXP2"
  ),
  
  Inhibitory_identity = c(
    "GAD1","GAD2","SLC32A1",
    "DLX1","DLX2","DLX5","DLX6",
    "LHX6","SOX6","ARX",
    "ERBB4","RELN"
  ),
  
  Inhibitory_subtypes = c(
    "SST","PVALB","VIP","CALB2","LAMP5"
  ),
  
  Presynaptic = c(
    "SYN1","SYN2","SYN3","SYP","SYT1","SYT2",
    "SNAP25","VAMP2","STX1A","STXBP1",
    "RIMS1","RIMS2","BSN","PCLO",
    "UNC13A","UNC13B","DNM1","DNM3"
  ),
  
  Postsynaptic = c(
    "DLG4","DLG3","DLG2",
    "SHANK1","SHANK2","SHANK3",
    "HOMER1","HOMER2","HOMER3",
    "SYNGAP1","ARC","NPTX1","NPTX2"
  ),
  
  Excitatory_synaptic_maturation = c(
    "CAMK2A","CAMK2B",
    "GRIA1","GRIA2","GRIA3",
    "GRIN1","GRIN2A","GRIN2B",
    "DLG4","SHANK1","SHANK3",
    "HOMER1","NPTX1","NPTX2"
  ),
  
  GABAergic_synaptic_maturation = c(
    "GAD1","GAD2","SLC32A1",
    "GABRA1","GABRA2","GABRA5",
    "GABRB1","GABRB2","GABRB3",
    "GABRG2","GPHN","ARHGEF9"
  ),
  
  Synaptic_vesicle_cycle = c(
    "SLC17A7","SLC17A6","SLC32A1",
    "SYN1","SYN2","SYP","SYT1",
    "VAMP2","SNAP25","STX1A","STXBP1",
    "RAB3A","RAB3B","RAB3C",
    "RIMS1","DNM1",
    "ATP6V1A","ATP6V1B2","ATP6V0C"
  ),
  
  Synaptic_adhesion = c(
    "NRXN1","NRXN2","NRXN3",
    "NLGN1","NLGN2","NLGN3",
    "CNTN1","CNTN2","CNTNAP2",
    "NCAM1","L1CAM","CHL1",
    "PTPRD","PTPRS"
  )
)
# ============================================================
# 4C.1 CALCULATE SYNAPSE MODULE SCORES
# ============================================================

synapse_scores <- purrr::map_dfc(names(synapse_modules), function(m) {
  tibble::tibble(!!m := score_module(mat_z, synapse_modules[[m]]))
}) %>%
  dplyr::mutate(sample = colnames(mat_z)) %>%
  dplyr::left_join(meta, by = "sample") %>%
  dplyr::relocate(sample, genotype, treatment, group)

write.xlsx(
  synapse_scores,
  file.path(output_dir, "Tables", "synapse_module_scores.xlsx")
)

synapse_detection <- purrr::imap_dfr(synapse_modules, function(genes, module) {
  tibble::tibble(
    module = module,
    gene = genes,
    detected = gene %in% rownames(mat_z)
  )
})

write.xlsx(
  synapse_detection,
  file.path(output_dir, "Tables", "synapse_module_gene_detection.xlsx")
)
# ============================================================
# 4D. SYNAPSE MODULE BOXPLOTS
# ============================================================

synapse_long <- synapse_scores %>%
  tidyr::pivot_longer(
    cols = names(synapse_modules),
    names_to = "Module",
    values_to = "Score"
  )

synapse_long$Module <- factor(synapse_long$Module, levels = names(synapse_modules))
synapse_long$group <- factor(
  synapse_long$group,
  levels = c(
    "NDC2_DMSO", "NDC2_DAPT",
    "sAD1_DMSO", "sAD1_DAPT",
    "sAD2_DMSO", "sAD2_DAPT"
  )
)

p_synapse_box <- ggplot(
  synapse_long,
  aes(x = group, y = Score, fill = treatment)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.80, width = 0.65) +
  geom_jitter(
    aes(shape = treatment),
    width = 0.10,
    size = 2.8,
    alpha = 0.9
  ) +
  facet_wrap(~Module, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 13) +
  labs(
    title = "Synaptic protein module scores",
    x = NULL,
    y = "Mean protein z-score"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    legend.position = "top"
  )

ggsave(
  file.path(output_dir, "Figures", "07_synapse_module_scores_boxplots.pdf"),
  p_synapse_box,
  width = 14,
  height = 7
)

ggsave(
  file.path(output_dir, "Figures", "07_synapse_module_scores_boxplots.png"),
  p_synapse_box,
  width = 14,
  height = 7,
  dpi = 300
)

# ============================================================
# 4E. SYNAPSE HEATMAP
# ============================================================

synapse_gene_map <- purrr::imap_dfr(synapse_modules, function(genes, module) {
  data.frame(Gene = genes, Module = module)
}) %>%
  dplyr::filter(Gene %in% rownames(mat_z)) %>%
  dplyr::mutate(Module = factor(Module, levels = names(synapse_modules))) %>%
  dplyr::arrange(Module, Gene) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

synapse_mat <- mat_z[synapse_gene_map$Gene, , drop = FALSE]

annotation_row_syn <- synapse_gene_map %>%
  dplyr::select(Gene, Module) %>%
  tibble::column_to_rownames("Gene")

annotation_col_syn <- meta %>%
  dplyr::select(sample, genotype, treatment) %>%
  tibble::column_to_rownames("sample")

synapse_module_colors <- list(
  Module = c(
    "Excitatory_identity" = "#E69F00",
    "Inhibitory_identity" = "#984EA3",
    "Inhibitory_subtypes" = "#7570B3",
    "Presynaptic" = "#1B7837",
    "Postsynaptic" = "#762A83",
    "Excitatory_synaptic_maturation" = "#D95F02",
    "GABAergic_synaptic_maturation" = "#4A7C59",
    "Synaptic_vesicle_cycle" = "#66A61E",
    "Synaptic_adhesion" = "#377EB8"
  ),
  
  treatment = c(
    "DMSO" = "#80B1D3",
    "DAPT" = "#C13C91"
  ),
  
  genotype = c(
    "NDC2" = "#E69F00",
    "sAD1" = "#984EA3",
    "sAD2" = "#006837"
  )
)
gaps_syn <- cumsum(table(synapse_gene_map$Module))
gaps_syn <- gaps_syn[-length(gaps_syn)]

pdf(
  file.path(output_dir, "Figures", "08_synapse_marker_heatmap_by_module.pdf"),
  width = 10,
  height = max(8, nrow(synapse_mat) * 0.22)
)

pheatmap(
  synapse_mat,
  annotation_col = annotation_col_syn,
  annotation_row = annotation_row_syn,
  annotation_colors = synapse_module_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_syn,
  color = colorRampPalette(c("#1B7837", "white", "#762A83"))(100),
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  main = "Proteomic synaptic markers"
)

dev.off()

# ============================================================
# 4F. CORRELATION BY LINE AND ENVIRONMENT:
# NEURON TYPE x SYNAPSE MODULES
# ============================================================

neuron_modules_for_corr <- c(
  "Excitatory_neurons",
  "Inhibitory_interneurons",
  "Intermediate_progenitors"
)

neuron_modules_for_corr <- intersect(neuron_modules_for_corr, colnames(module_scores))

corr_input <- module_scores %>%
  dplyr::select(sample, genotype, treatment, group, dplyr::all_of(neuron_modules_for_corr)) %>%
  dplyr::left_join(
    synapse_scores %>%
      dplyr::select(sample, dplyr::all_of(names(synapse_modules))),
    by = "sample"
  )

corr_results_by_group <- corr_input %>%
  dplyr::group_by(genotype, treatment, group) %>%
  dplyr::group_modify(~{
    
    df <- .x
    
    expand.grid(
      Neuronal_module = neuron_modules_for_corr,
      Synapse_module = names(synapse_modules),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n = sum(
          complete.cases(
            df[[Neuronal_module]],
            df[[Synapse_module]]
          )
        ),
        cor = ifelse(
          n >= 3,
          cor(
            df[[Neuronal_module]],
            df[[Synapse_module]],
            method = "spearman",
            use = "pairwise.complete.obs"
          ),
          NA_real_
        ),
        p_value = ifelse(
          n >= 3,
          suppressWarnings(
            cor.test(
              df[[Neuronal_module]],
              df[[Synapse_module]],
              method = "spearman"
            )$p.value
          ),
          NA_real_
        )
      ) %>%
      dplyr::ungroup()
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()

write.xlsx(
  corr_results_by_group,
  file.path(output_dir, "Tables", "neuronal_type_vs_synapse_correlations_by_line_treatment.xlsx")
)

# ============================================================
# 4F.1 HEATMAP — ONE PANEL PER LINE/ENVIRONMENT
# ============================================================

corr_results_by_group$Facet <- paste(
  corr_results_by_group$genotype,
  corr_results_by_group$treatment,
  sep = " - "
)

p_corr_by_group <- ggplot(
  corr_results_by_group,
  aes(
    x = Synapse_module,
    y = Neuronal_module,
    fill = cor
  )
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(cor), "", sprintf("%.2f", cor))),
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(~Facet, ncol = 3) +
  scale_fill_gradient2(
    low = "#1B7837",
    mid = "white",
    high = "#762A83",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Neuron subtype–synapse correlations by line and treatment",
    x = NULL,
    y = NULL,
    fill = "Spearman r"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "Figures", "09_neuronal_type_vs_synapse_correlation_by_line_treatment.pdf"),
  p_corr_by_group,
  width = 13,
  height = 7
)

ggsave(
  file.path(output_dir, "Figures", "09_neuronal_type_vs_synapse_correlation_by_line_treatment.png"),
  p_corr_by_group,
  width = 13,
  height = 7,
  dpi = 300
)
  # ============================================================
# CLEAN CORRELATION: NEURONAL IDENTITY x SYNAPTIC FUNCTION
# ============================================================

neuron_modules_for_corr <- c(
  "Excitatory_identity",
  "Inhibitory_identity"
)

synapse_modules_for_corr <- c(
  "Presynaptic",
  "Postsynaptic",
  "Synaptic_vesicle_cycle",
  "Excitatory_synaptic_maturation",
  "GABAergic_synaptic_maturation"
)

corr_input_clean <- synapse_scores %>%
  dplyr::select(
    sample, genotype, treatment, group,
    dplyr::all_of(neuron_modules_for_corr),
    dplyr::all_of(synapse_modules_for_corr)
  )

corr_clean_by_group <- corr_input_clean %>%
  dplyr::group_by(genotype, treatment, group) %>%
  dplyr::group_modify(~{

    df <- .x

    expand.grid(
      Neuronal_module = neuron_modules_for_corr,
      Synapse_module = synapse_modules_for_corr,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n = sum(complete.cases(df[[Neuronal_module]], df[[Synapse_module]])),
        cor = ifelse(
          n >= 3,
          cor(
            df[[Neuronal_module]],
            df[[Synapse_module]],
            method = "spearman",
            use = "pairwise.complete.obs"
          ),
          NA_real_
        )
      ) %>%
      dplyr::ungroup()
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Facet = paste(genotype, treatment, sep = " - ")
  )

write.xlsx(
  corr_clean_by_group,
  file.path(output_dir, "Tables", "clean_neuronal_identity_synapse_correlations_by_group.xlsx")
)

p_corr_clean <- ggplot(
  corr_clean_by_group,
  aes(x = Synapse_module, y = Neuronal_module, fill = cor)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = ifelse(is.na(cor), "", sprintf("%.2f", cor))),
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~Facet, ncol = 3) +
  scale_fill_gradient2(
    low = "#1B7837",
    mid = "white",
    high = "#762A83",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Correlation between neuronal identity and synaptic maturation",
    x = NULL,
    y = NULL,
    fill = "Spearman r"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "Figures", "09_clean_identity_vs_synapse_correlation_by_group.pdf"),
  p_corr_clean,
  width = 11,
  height = 5.5
)

ggsave(
  file.path(output_dir, "Figures", "09_clean_identity_vs_synapse_correlation_by_group.png"),
  p_corr_clean,
  width = 11,
  height = 5.5,
  dpi = 300
)
# ============================================================
# CLEAN CORRELATION: NEURONAL IDENTITY x SYNAPTIC FUNCTION
# ============================================================

neuron_modules_for_corr <- c(
  "Excitatory_identity",
  "Inhibitory_identity"
)

synapse_modules_for_corr <- c(
  "Presynaptic",
  "Postsynaptic",
  "Synaptic_vesicle_cycle",
  "Excitatory_synaptic_maturation",
  "GABAergic_synaptic_maturation"
)

corr_input_clean <- synapse_scores %>%
  dplyr::select(
    sample, genotype, treatment, group,
    dplyr::all_of(neuron_modules_for_corr),
    dplyr::all_of(synapse_modules_for_corr)
  )

corr_clean_by_group <- corr_input_clean %>%
  dplyr::group_by(genotype, treatment, group) %>%
  dplyr::group_modify(~{
    
    df <- .x
    
    expand.grid(
      Neuronal_module = neuron_modules_for_corr,
      Synapse_module = synapse_modules_for_corr,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n = sum(complete.cases(df[[Neuronal_module]], df[[Synapse_module]])),
        cor = ifelse(
          n >= 3,
          cor(
            df[[Neuronal_module]],
            df[[Synapse_module]],
            method = "spearman",
            use = "pairwise.complete.obs"
          ),
          NA_real_
        )
      ) %>%
      dplyr::ungroup()
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Facet = paste(genotype, treatment, sep = " - ")
  )

write.xlsx(
  corr_clean_by_group,
  file.path(output_dir, "Tables", "clean_neuronal_identity_synapse_correlations_by_group.xlsx")
)

p_corr_clean <- ggplot(
  corr_clean_by_group,
  aes(x = Synapse_module, y = Neuronal_module, fill = cor)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = ifelse(is.na(cor), "", sprintf("%.2f", cor))),
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(~Facet, ncol = 3) +
  scale_fill_gradient2(
    low = "#1B7837",
    mid = "white",
    high = "#762A83",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Correlation between neuronal identity and synaptic maturation",
    x = NULL,
    y = NULL,
    fill = "Spearman r"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(output_dir, "Figures", "09_clean_identity_vs_synapse_correlation_by_group.pdf"),
  p_corr_clean,
  width = 11,
  height = 5.5
)

ggsave(
  file.path(output_dir, "Figures", "09_clean_identity_vs_synapse_correlation_by_group.png"),
  p_corr_clean,
  width = 11,
  height = 5.5,
  dpi = 300
)
# ============================================================
# DELTA CORRELATION (DAPT - DMSO)
# ============================================================

delta_corr <- corr_clean_by_group %>%
  
  dplyr::select(
    genotype,
    treatment,
    Neuronal_module,
    Synapse_module,
    cor
  ) %>%
  
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = cor
  ) %>%
  
  dplyr::mutate(
    Delta = DAPT - DMSO
  )

delta_heatmap <- delta_corr %>%
  
  dplyr::mutate(
    Pair = paste(Neuronal_module, Synapse_module, sep = " → ")
  ) %>%
  
  dplyr::select(
    genotype,
    Pair,
    Delta
  ) %>%
  
  tidyr::pivot_wider(
    names_from = genotype,
    values_from = Delta
  ) %>%
  
  tibble::column_to_rownames("Pair") %>%
  
  as.matrix()

pdf(
  file.path(output_dir,
  "Figures",
  "10_delta_correlation_DAPT_minus_DMSO.pdf"),
  width = 5.5,
  height = 7
)

pheatmap(
  delta_heatmap,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  color = colorRampPalette(
    c("#1B7837", "white", "#762A83")
  )(100),
  
  breaks = seq(-2, 2, length.out = 101),
  
  fontsize_row = 9,
  fontsize_col = 11,
  
  border_color = NA,
  
  display_numbers = TRUE,
  number_format = "%.2f",
  
  main = "Δ correlation (DAPT − DMSO)"
)
dev.off()


# ============================================================
# SYNAPTIC MATURATION INDEX
# ============================================================

synaptic_maturation_index <- synapse_scores %>%
  dplyr::mutate(
    Synaptic_Maturation_Index =
      rowMeans(
        dplyr::select(
          .,
          dplyr::any_of(c(
            "Presynaptic",
            "Postsynaptic",
            "Excitatory_synaptic_maturation",
            "GABAergic_synaptic_maturation",
            "Synaptic_vesicle_cycle",
            "Synaptic_adhesion"
          ))
        ),
        na.rm = TRUE
      ),
    
    E_I_Balance_Index =
      Excitatory_synaptic_maturation -
      GABAergic_synaptic_maturation
  )
# =========================
# 5. BOXPLOTS OF MODULE SCORES
# =========================

module_scores$group <- factor(
  module_scores$group,
  levels = c(
    "NDC2_DMSO", "NDC2_DAPT",
    "sAD1_DMSO", "sAD1_DAPT",
    "sAD2_DMSO", "sAD2_DAPT"
  )
)

score_long <- module_scores %>%
  tidyr::pivot_longer(
    cols = names(modules),
    names_to = "Module",
    values_to = "Score"
  )

score_long$Module <- factor(
  score_long$Module,
  levels = names(modules)
)

group_colors <- c(
  "DMSO" = "#80B1D3",
  "DAPT" = "#C13C91"
)

p_modules <- ggplot(
  score_long,
  aes(x = group, y = Score, fill = treatment)
) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.80,
    width = 0.65
  ) +
  
  geom_jitter(
    aes(shape = treatment),
    width = 0.10,
    size = 2.8,
    alpha = 0.9
  ) +
  
  facet_wrap(
    ~Module,
    scales = "free_y",
    ncol = 3
  ) +
  
  scale_fill_manual(values = group_colors) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "Cell-type and neuronal maturation module scores",
    x = NULL,
    y = "Mean protein z-score"
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(face = "bold"),
    
    axis.title.y = element_text(face = "bold"),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    
    legend.position = "top"
  )

ggsave(
  file.path(output_dir, "Figures", "05_module_scores_boxplots.pdf"),
  p_modules,
  width = 14,
  height = 10
)

ggsave(
  file.path(output_dir, "Figures", "05_module_scores_boxplots.png"),
  p_modules,
  width = 14,
  height = 10,
  dpi = 300
)
# =========================
# 6. HEATMAP OF MODULE SCORES
# =========================

score_mat <- module_scores %>%
  dplyr::select(sample, dplyr::all_of(names(modules))) %>%
  tibble::column_to_rownames("sample") %>%
  as.matrix()

# Remove módulos que ficaram completamente NA
score_mat <- score_mat[, colSums(!is.na(score_mat)) > 0, drop = FALSE]

# Scale
score_mat_scaled <- scale(score_mat)

# Troca NA/NaN/Inf por 0 para o pheatmap não quebrar
score_mat_scaled[is.na(score_mat_scaled)] <- 0
score_mat_scaled[is.nan(score_mat_scaled)] <- 0
score_mat_scaled[is.infinite(score_mat_scaled)] <- 0

annotation_col <- meta %>%
  dplyr::select(sample, genotype, treatment, group) %>%
  tibble::column_to_rownames("sample")

pdf(file.path(output_dir, "Figures", "02_module_score_heatmap.pdf"),
    width = 8,
    height = 6)

pheatmap(
  t(score_mat_scaled),
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#1B7837", "white", "#762A83"))(100),
  fontsize = 10,
  main = "Proteomic neuronal maturation module scores"
)

dev.off()

# =========================
# 7. CORE PROTEIN HEATMAP BY CELL-TYPE MODULES
# =========================

module_order <- names(modules)

gene_module_map <- purrr::imap_dfr(modules, function(genes, module) {
  data.frame(Gene = genes, Module = module)
}) %>%
  dplyr::filter(Gene %in% rownames(mat_z)) %>%
  dplyr::mutate(Module = factor(Module, levels = module_order)) %>%
  dplyr::arrange(Module, Gene) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

z_mat_module <- mat_z[gene_module_map$Gene, , drop = FALSE]

annotation_row <- gene_module_map %>%
  dplyr::select(Gene, Module) %>%
  tibble::column_to_rownames("Gene")

annotation_col <- meta %>%
  dplyr::select(sample, genotype, treatment) %>%
  tibble::column_to_rownames("sample")

module_colors <- list(
  Module = c(
    "Radial_glia" = "#8B1E3F",
    "Cycling_progenitors" = "#C13C91",
    "Intermediate_progenitors" = "#D95F02",
    "Excitatory_neurons" = "#1B7837",
    "Inhibitory_interneurons" = "#4A7C59",
    "Astrocyte" = "#1F78B4",
    "Reactive_astrocyte" = "#762A83",
    "Oligodendrocyte_lineage" = "#377EB8",
    "Choroid_plexus" = "#E6AB02",
    "Mesenchymal_like" = "#A6761D"
  ),
  
  treatment = c(
    "DMSO" = "#80B1D3",
    "DAPT" = "#C13C91"
  ),
  
  genotype = c(
    "NDC2" = "#E69F00",
    "sAD1" = "#984EA3",
    "sAD2" = "#006837"
  )
)

gaps_row <- cumsum(table(gene_module_map$Module))
gaps_row <- gaps_row[-length(gaps_row)]

pdf(
  file.path(output_dir, "Figures", "03_core_celltype_maturation_heatmap_by_module_clean.pdf"),
  width = 10,
  height = max(8, nrow(z_mat_module) * 0.20)
)

pheatmap(
  z_mat_module,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = module_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_row = gaps_row,
  color = colorRampPalette(c("#1B7837", "white", "#762A83"))(100),
  fontsize_row = 7,
  fontsize_col = 8,
  border_color = NA,
  main = "Proteomic cell-type and neuronal maturation markers"
)

dev.off()

core_present <- unique(gene_module_map$Gene)
# =========================
# 8. PCA USING MATURATION PROTEINS
# =========================

pca <- prcomp(t(mat_z[core_present, , drop = FALSE]), center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::left_join(meta, by = "sample")

var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = genotype, shape = treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_classic(base_size = 13) +
  labs(
    title = "PCA based on neuronal differentiation/maturation proteins",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  )

ggsave(file.path(output_dir, "Figures", "04_PCA_neuronal_maturation_proteins.pdf"),
       p_pca, width = 6, height = 5)

ggsave(file.path(output_dir, "Figures", "04_PCA_neuronal_maturation_proteins.png"),
       p_pca, width = 6, height = 5, dpi = 300)

# =========================
# 9. LIMMA DIFFERENTIAL ANALYSIS
# =========================

meta$group_limma <- make.names(as.character(meta$group))

design <- model.matrix(~0 + group_limma, data = meta)
colnames(design) <- levels(factor(meta$group_limma))

fit <- lmFit(mat_imp, design)

contrast_matrix <- makeContrasts(
  NDC2_DAPT_vs_DMSO = NDC2_DAPT - NDC2_DMSO,
  sAD1_DAPT_vs_DMSO = sAD1_DAPT - sAD1_DMSO,
  sAD2_DAPT_vs_DMSO = sAD2_DAPT - sAD2_DMSO,
  
  sAD1_DMSO_vs_NDC2_DMSO = sAD1_DMSO - NDC2_DMSO,
  sAD2_DMSO_vs_NDC2_DMSO = sAD2_DMSO - NDC2_DMSO,
  
  sAD1_DAPT_vs_NDC2_DAPT = sAD1_DAPT - NDC2_DAPT,
  sAD2_DAPT_vs_NDC2_DAPT = sAD2_DAPT - NDC2_DAPT,
  
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

limma_results <- list()

for (cn in colnames(contrast_matrix)) {
  res <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH") %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::mutate(Contrast = cn)
  
  limma_results[[cn]] <- res
  
  write.xlsx(res, file.path(output_dir, "Tables", paste0("limma_", cn, ".xlsx")))
}

all_limma <- dplyr::bind_rows(limma_results)
write.xlsx(all_limma, file.path(output_dir, "Tables", "ALL_limma_results.xlsx"))


# =========================
# 11. DAPT MATURATION INDEX
# =========================

positive_modules <- c(
  "Intermediate_progenitors",
  "Excitatory_neurons",
  "Inhibitory_interneurons"
)

negative_modules <- c(
  "Radial_glia",
  "Cycling_progenitors",
  "Reactive_astrocyte"
)

positive_modules <- intersect(positive_modules, colnames(module_scores))
negative_modules <- intersect(negative_modules, colnames(module_scores))

maturation_index <- module_scores %>%
  dplyr::mutate(
    Maturation_Index =
      rowMeans(dplyr::select(., dplyr::all_of(positive_modules)), na.rm = TRUE) -
      rowMeans(dplyr::select(., dplyr::all_of(negative_modules)), na.rm = TRUE)
  )

write.xlsx(
  maturation_index,
  file.path(output_dir, "Tables", "maturation_index.xlsx")
)

p_index <- ggplot(maturation_index, aes(x = group, y = Maturation_Index, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
  geom_jitter(aes(shape = treatment), width = 0.12, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = c("DMSO" = "#80B1D3", "DAPT" = "#C13C91")) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = NULL,
    y = "Proteomic maturation index",
    title = "Proteomic maturation index"
  )

ggsave(file.path(output_dir, "Figures", "06_maturation_index.pdf"),
       p_index, width = 7, height = 5)

ggsave(file.path(output_dir, "Figures", "06_maturation_index.png"),
       p_index, width = 7, height = 5, dpi = 300)
# =========================
# 13. ADDITIONAL MATURATION FIGURES
# =========================

extra_fig_dir <- file.path(output_dir, "Figures", "Additional_maturation_plots")
dir.create(extra_fig_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(
  "DMSO" = "#80B1D3",
  "DAPT" = "#C13C91"
)

# -------------------------
# A. MODULE TRAJECTORY PLOT
# -------------------------

trajectory_modules <- c(
  "Radial_glia",
  "Cycling_progenitors",
  "Intermediate_progenitors",
  "Excitatory_neurons",
  "Inhibitory_interneurons",
  "Oligodendrocyte_lineage"
)

trajectory_modules <- intersect(trajectory_modules, colnames(module_scores))

trajectory_df <- module_scores %>%
  dplyr::select(sample, genotype, treatment, group, dplyr::all_of(trajectory_modules)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(trajectory_modules),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  dplyr::group_by(genotype, treatment, group, Module) %>%
  dplyr::summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Module = factor(Module, levels = trajectory_modules)
  )

p_traj <- ggplot(
  trajectory_df,
  aes(x = Module, y = mean_score, group = treatment, color = treatment)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
    width = 0.15,
    linewidth = 0.6
  ) +
  facet_wrap(~genotype, nrow = 1) +
  scale_color_manual(values = group_colors) +
  theme_classic(base_size = 13) +
  labs(
    title = "Proteomic developmental trajectory across cell-state modules",
    x = NULL,
    y = "Mean module z-score"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(extra_fig_dir, "A_module_trajectory_plot.pdf"),
       p_traj, width = 11, height = 5)

ggsave(file.path(extra_fig_dir, "A_module_trajectory_plot.png"),
       p_traj, width = 11, height = 5, dpi = 300)


# -------------------------
# B. PROGENITOR VS MATURATION SCATTER
# -------------------------

scatter_df <- module_scores %>%
  dplyr::mutate(
    Progenitor_score = rowMeans(
      dplyr::select(
        .,
        dplyr::any_of(c("Radial_glia", "Cycling_progenitors"))
      ),
      na.rm = TRUE
    ),
    Maturation_score = rowMeans(
      dplyr::select(
        .,
        dplyr::any_of(c(
          "Intermediate_progenitors",
          "Excitatory_neurons",
          "Inhibitory_interneurons",
          "Oligodendrocyte_lineage"
        ))
      ),
      na.rm = TRUE
    )
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = Progenitor_score,
    y = Maturation_score,
    color = treatment,
    shape = genotype
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(
    aes(group = interaction(genotype, treatment)),
    linewidth = 0.5,
    alpha = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  theme_classic(base_size = 13) +
  labs(
    title = "Proteomic shift from progenitor to maturation state",
    x = "Progenitor module score",
    y = "Maturation module score"
  ) +
  theme(
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(extra_fig_dir, "B_progenitor_vs_maturation_scatter.pdf"),
       p_scatter, width = 6.5, height = 5.5)

ggsave(file.path(extra_fig_dir, "B_progenitor_vs_maturation_scatter.png"),
       p_scatter, width = 6.5, height = 5.5, dpi = 300)


# -------------------------
# C. DELTA PLOT: DAPT - GROUND
# -------------------------

delta_modules <- setdiff(
  names(modules),
  c("Endothelial", "Caudal_neurons")
)

delta_modules <- intersect(delta_modules, colnames(module_scores))

delta_df <- module_scores %>%
  dplyr::select(sample, genotype, treatment, dplyr::all_of(delta_modules)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(delta_modules),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  dplyr::group_by(genotype, treatment, Module) %>%
  dplyr::summarise(
    mean_score = mean(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = treatment,
    values_from = mean_score
  ) %>%
  dplyr::mutate(
    Delta_DAPT_minus_DMSO = DAPT - DMSO,
    Module = factor(Module, levels = delta_modules)
  )

write.xlsx(
  delta_df,
  file.path(output_dir, "Tables", "delta_DAPT_minus_DMSO_by_module.xlsx")
)

p_delta <- ggplot(
  delta_df,
  aes(x = Module, y = Delta_DAPT_minus_DMSO, fill = genotype)
) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7,
    color = "black",
    linewidth = 0.3,
    alpha = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "NDC2" = "#E69F00",
      "sAD1" = "#984EA3",
      "sAD2" = "#006837"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Spaceflight-associated module shift",
    x = NULL,
    y = "Δ module score (DAPT - DMSO)"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(extra_fig_dir, "C_delta_DAPT_minus_DMSO_by_module.pdf"),
       p_delta, width = 11, height = 5.5)

ggsave(file.path(extra_fig_dir, "C_delta_DAPT_minus_DMSO_by_module.png"),
       p_delta, width = 11, height = 5.5, dpi = 300)


# -------------------------
# D. DELTA HEATMAP
# -------------------------

delta_mat <- delta_df %>%
  dplyr::select(genotype, Module, Delta_DAPT_minus_DMSO) %>%
  tidyr::pivot_wider(
    names_from = genotype,
    values_from = Delta_DAPT_minus_DMSO
  ) %>%
  tibble::column_to_rownames("Module") %>%
  as.matrix()

pdf(file.path(extra_fig_dir, "D_delta_DAPT_minus_DMSO_heatmap.pdf"),
    width = 5.5,
    height = max(5, nrow(delta_mat) * 0.35))

pheatmap(
  delta_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#1B7837", "white", "#762A83"))(100),
  fontsize = 10,
  border_color = NA,
  main = "Δ module score: DAPT - DMSO"
)

dev.off()
# =========================
# 12. EXPORT SUMMARY OF DETECTED MARKERS
# =========================

summary_detected <- module_detection %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_detected = sum(detected),
    detected_genes = paste(gene[detected], collapse = ", "),
    missing_genes = paste(gene[!detected], collapse = ", "),
    .groups = "drop"
  )

write.xlsx(summary_detected,
           file.path(output_dir, "Tables", "module_detection_summary.xlsx"))

cat("\nAnalysis completed successfully!\n")
cat("All files saved in:", normalizePath(output_dir), "\n")