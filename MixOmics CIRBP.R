install.packages("BiocManager")
BiocManager::install("mixOmics")
BiocManager::install("mixOmicsTeam/mixOmics")
install.packages("ppcor")
install.packages("mediation")

###########################################################################

## =========================
## 0) PACKAGES (install if missing)
## =========================
pkgs_cran <- c("openxlsx","pheatmap","ggplot2","emmeans","dplyr")
pkgs_bioc <- c("mixOmics","BiocParallel","limma")

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
}
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
have_ragg <- requireNamespace("ragg", quietly=TRUE)

suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))

## recommended stability
register(SerialParam())

## =========================
## 1) INPUTS / OUTPUTS
## =========================
rna_path  <- "/Users/namoreira/Downloads/Transcriptomic spx30.csv"
prot_path <- "/Users/namoreira/Downloads/PROTEOMICS SPX30.csv"
outdir    <- "/Users/namoreira/Downloads/mixomics_DIABLO_PLANB_CLEAN_UPDATED"

dir_fig_pdf  <- file.path(outdir, "figures_pdf")
dir_fig_tiff <- file.path(outdir, "figures_tiff")
dir_tables   <- file.path(outdir, "tables_csv")

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(dir_fig_pdf, showWarnings=FALSE, recursive=TRUE)
dir.create(dir_fig_tiff, showWarnings=FALSE, recursive=TRUE)
dir.create(dir_tables, showWarnings=FALSE, recursive=TRUE)

cat("Output folder:\n", outdir, "\n\n", sep="")

## =========================
## 2) SAFE HELPERS
## =========================
safe_dev_off <- function() {
  while (grDevices::dev.cur() > 1) {
    try(grDevices::dev.off(), silent=TRUE)
  }
}

save_plot_pdf <- function(filename, plot_expr, width=10, height=8) {
  safe_dev_off()
  outfile <- file.path(dir_fig_pdf, filename)
  ok <- TRUE
  grDevices::pdf(outfile, width=width, height=height, useDingbats=FALSE)
  tryCatch({ force(plot_expr) },
           error=function(e){ ok <<- FALSE; message("\n[PDF FAILED] ", filename, "\n", conditionMessage(e), "\n") },
           finally={ try(grDevices::dev.off(), silent=TRUE) })
  invisible(ok)
}

save_plot_tiff <- function(filename, plot_expr, width_in=10, height_in=8, dpi=600) {
  safe_dev_off()
  outfile <- file.path(dir_fig_tiff, filename)
  w_px <- as.integer(width_in * dpi)
  h_px <- as.integer(height_in * dpi)
  ok <- TRUE
  tryCatch({
    if (have_ragg) {
      ragg::agg_tiff(outfile, width=w_px, height=h_px, res=dpi)
    } else {
      grDevices::tiff(outfile, width=w_px, height=h_px, res=dpi, compression="lzw")
    }
    force(plot_expr)
  }, error=function(e){
    ok <<- FALSE
    message("\n[TIFF FAILED] ", filename, "\n", conditionMessage(e), "\n")
  }, finally={
    try(grDevices::dev.off(), silent=TRUE)
  })
  invisible(ok)
}

write_csv <- function(df, filename) {
  outfile <- file.path(dir_tables, filename)
  utils::write.csv(df, outfile, row.names=FALSE)
  invisible(outfile)
}

sanitize_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

## =========================
## 2B) READ + ROBUST NUMERIC CONVERSION
## Expectation: first column = feature IDs, remaining columns = samples
## =========================
read_matrix_csv <- function(path) {
  df <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  if (ncol(df) < 2) stop("File has <2 columns: ", path)
  feat <- df[[1]]
  df  <- df[, -1, drop=FALSE]
  
  # Convert all to numeric safely (handles "1,23" if needed)
  to_num <- function(x) {
    if (is.numeric(x)) return(x)
    x <- gsub(",", ".", x)
    suppressWarnings(as.numeric(x))
  }
  for (j in seq_len(ncol(df))) df[[j]] <- to_num(df[[j]])
  
  mat <- as.matrix(df)
  rownames(mat) <- feat
  mat
}

## =========================
## 2C) PREPROCESS HELPERS
## =========================
impute_median <- function(x) {
  for (j in seq_len(ncol(x))) {
    nas <- is.na(x[, j])
    if (any(nas)) x[nas, j] <- stats::median(x[, j], na.rm=TRUE)
  }
  x
}

as_nzv_index <- function(nzv_out, p) {
  if (is.null(nzv_out)) return(integer(0))
  if (is.data.frame(nzv_out)) {
    if ("nzv" %in% names(nzv_out)) nzv_out <- which(nzv_out$nzv %in% TRUE)
    else if ("nearZeroVar" %in% names(nzv_out)) nzv_out <- which(nzv_out$nearZeroVar %in% TRUE)
    else if ("zeroVar" %in% names(nzv_out)) nzv_out <- which(nzv_out$zeroVar %in% TRUE)
    else return(integer(0))
  }
  if (is.logical(nzv_out)) {
    if (length(nzv_out) != p) return(integer(0))
    nzv_out <- which(nzv_out %in% TRUE)
  }
  if (is.list(nzv_out) && !is.numeric(nzv_out)) {
    if (!is.null(nzv_out$nzv)) nzv_out <- nzv_out$nzv
    else if (!is.null(nzv_out$index)) nzv_out <- nzv_out$index
    else return(integer(0))
  }
  if (!is.numeric(nzv_out)) return(integer(0))
  nzv_out <- as.integer(nzv_out)
  nzv_out <- nzv_out[!is.na(nzv_out)]
  nzv_out <- nzv_out[nzv_out >= 1 & nzv_out <= p]
  unique(nzv_out)
}

remove_nzv <- function(x) {
  p <- ncol(x)
  if (p == 0) return(x)
  nzv_raw <- try(mixOmics::nearZeroVar(x), silent=TRUE)
  if (inherits(nzv_raw, "try-error")) return(x)
  idx <- as_nzv_index(nzv_raw, p)
  if (length(idx) == 0) return(x)
  keep <- setdiff(seq_len(p), idx)
  if (length(keep) == 0) return(x[, 0, drop=FALSE])
  x[, keep, drop=FALSE]
}

## =========================
## 2D) HEATMAP (NO sample clustering; keeps sample order)
## group must be factor length == nrow(mat)
## =========================
make_heatmap_nocluster <- function(mat, group, filename_prefix,
                                   scale_rows=TRUE,
                                   show_rownames=FALSE,
                                   fontsize_col=10,
                                   cluster_features=TRUE) {
  
  if (is.null(mat) || nrow(mat) < 2 || ncol(mat) < 1) {
    cat("Skipping heatmap ", filename_prefix, " (dim=", nrow(mat), "x", ncol(mat), ")\n", sep="")
    return(invisible(NULL))
  }
  stopifnot(nrow(mat) == length(group))
  group <- factor(group)
  
  if (ncol(mat) < 2) cluster_features <- FALSE
  
  mat_plot <- mat
  if (scale_rows) {
    mat_plot <- t(scale(t(mat_plot)))
    mat_plot[is.na(mat_plot)] <- 0
  }
  
  ann <- data.frame(Group=group)
  rownames(ann) <- rownames(mat_plot)
  
  grp_cols <- setNames(mixOmics::color.mixo(seq_along(levels(group))), levels(group))
  ann_cols <- list(Group=grp_cols)
  
  save_plot_pdf(paste0(filename_prefix, ".pdf"), {
    ph <- pheatmap::pheatmap(
      mat_plot,
      cluster_rows=FALSE,
      cluster_cols=cluster_features,
      annotation_row=ann,
      annotation_colors=ann_cols,
      show_rownames=show_rownames,
      fontsize_col=fontsize_col,
      border_color=NA,
      silent=TRUE
    )
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
  }, width=12, height=8)
  
  save_plot_tiff(paste0(filename_prefix, ".tiff"), {
    ph <- pheatmap::pheatmap(
      mat_plot,
      cluster_rows=FALSE,
      cluster_cols=cluster_features,
      annotation_row=ann,
      annotation_colors=ann_cols,
      show_rownames=show_rownames,
      fontsize_col=fontsize_col,
      border_color=NA,
      silent=TRUE
    )
    grid::grid.newpage(); grid::grid.draw(ph$gtable)
  }, width_in=12, height_in=8, dpi=600)
  
  invisible(TRUE)
}

## =========================
## 2E) mixOmics selectVar extraction (robust)
## =========================
get_selected_features <- function(model, block, comp) {
  flatten_chr <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.list(x)) x <- unlist(x, use.names=FALSE)
    as.character(x)
  }
  flatten_num <- function(x) {
    if (is.null(x)) return(numeric(0))
    if (is.list(x)) x <- unlist(x, use.names=FALSE)
    suppressWarnings(as.numeric(x))
  }
  
  out <- try(mixOmics::selectVar(model, block=block, comp=comp), silent=TRUE)
  if (!inherits(out, "try-error") && is.list(out) && (!is.null(out$name) || !is.null(out$value))) {
    nm  <- flatten_chr(out$name)
    val <- flatten_num(out$value)
    L <- min(length(nm), length(val))
    if (L == 0) return(data.frame(block=character(0), comp=integer(0), feature=character(0), loading=numeric(0)))
    return(data.frame(block=rep(block,L), comp=rep(comp,L), feature=nm[1:L], loading=val[1:L], stringsAsFactors=FALSE))
  }
  
  out2 <- try(mixOmics::selectVar(model, comp=comp), silent=TRUE)
  if (inherits(out2, "try-error") || !is.list(out2) || is.null(out2[[block]])) {
    return(data.frame(block=character(0), comp=integer(0), feature=character(0), loading=numeric(0)))
  }
  b <- out2[[block]]
  nm <- flatten_chr(b$name); val <- flatten_num(b$value)
  L <- min(length(nm), length(val))
  if (L == 0) return(data.frame(block=character(0), comp=integer(0), feature=character(0), loading=numeric(0)))
  data.frame(block=rep(block,L), comp=rep(comp,L), feature=nm[1:L], loading=val[1:L], stringsAsFactors=FALSE)
}

make_cross_cor <- function(model, X, comp=1, block1="transcriptome", block2="proteome") {
  sel1 <- try(mixOmics::selectVar(model, block=block1, comp=comp)$name, silent=TRUE)
  sel2 <- try(mixOmics::selectVar(model, block=block2, comp=comp)$name, silent=TRUE)
  if (inherits(sel1,"try-error") || inherits(sel2,"try-error")) {
    return(list(matrix=matrix(numeric(0),0,0), long=data.frame()))
  }
  sel1 <- sel1[!is.na(sel1)]
  sel2 <- sel2[!is.na(sel2)]
  if (length(sel1)==0 || length(sel2)==0) return(list(matrix=matrix(numeric(0),0,0), long=data.frame()))
  
  m1 <- X[[block1]][, sel1, drop=FALSE]
  m2 <- X[[block2]][, sel2, drop=FALSE]
  cmat <- stats::cor(m1, m2, use="pairwise.complete.obs", method="pearson")
  
  long <- data.frame(
    feature_block1 = rep(rownames(cmat), times=ncol(cmat)),
    feature_block2 = rep(colnames(cmat), each=nrow(cmat)),
    correlation    = as.vector(cmat),
    stringsAsFactors=FALSE
  )
  list(matrix=cmat, long=long)
}

## =========================
## 3) LOAD + TRANSPOSE to samples x features
## =========================
rna_raw  <- read_matrix_csv(rna_path)
prot_raw <- read_matrix_csv(prot_path)

rna  <- t(rna_raw)   # samples x genes
prot <- t(prot_raw)  # samples x proteins

cat("RNA dim  (samples x genes): ", paste(dim(rna), collapse=" x "), "\n", sep="")
cat("PROT dim (samples x prots): ", paste(dim(prot), collapse=" x "), "\n\n", sep="")

## =========================
## 4) ALIGN SAMPLES (keep RNA order)
## =========================
common <- rownames(rna)[rownames(rna) %in% rownames(prot)]
if (length(common) < 6) stop("Too few common samples. Check sample names between RNA and Proteomics files.")

rna  <- rna[common, , drop=FALSE]
prot <- prot[common, , drop=FALSE]
stopifnot(identical(rownames(rna), rownames(prot)))

## =========================
## 5) PREPROCESS
## =========================
# drop zero-variance
rna  <- rna[,  apply(rna,  2, stats::var, na.rm=TRUE) > 0, drop=FALSE]
prot <- prot[, apply(prot, 2, stats::var, na.rm=TRUE) > 0, drop=FALSE]

# impute
rna  <- impute_median(rna)
prot <- impute_median(prot)

# scale
rna  <- scale(rna)
prot <- scale(prot)

# remove NZV
rna  <- remove_nzv(rna)
prot <- remove_nzv(prot)

## =========================
## 6) DEFINE META (from sample names)
## =========================
s <- rownames(rna)

geno <- ifelse(grepl("WT", s, ignore.case=TRUE), "WT",
               ifelse(grepl("APP",  s, ignore.case=TRUE), "APP",
                      ifelse(grepl("PSEN1",s, ignore.case=TRUE), "PSEN1",
                             ifelse(grepl("PSEN2",s, ignore.case=TRUE), "PSEN2", NA))))

env <- ifelse(grepl("GROUND", s, ignore.case=TRUE), "GROUND",
              ifelse(grepl("ISS",    s, ignore.case=TRUE), "ISS", NA))

geno <- factor(geno)
env  <- factor(env, levels=c("GROUND","ISS"))

Y     <- factor(paste(geno, env, sep="_"))
Y_env <- env

if (any(is.na(Y))) {
  bad <- s[is.na(Y)]
  stop("Some samples could not be assigned to genotype/environment. Examples:\n",
       paste(head(bad, 10), collapse=", "))
}

cat("Group counts (Genotype_Environment):\n"); print(table(Y)); cat("\n")
cat("Group counts (Environment only):\n");      print(table(Y_env)); cat("\n")

meta <- data.frame(
  sample = s,
  genotype = geno,
  environment = env,
  group = Y,
  stringsAsFactors=FALSE
)
rownames(meta) <- meta$sample

## =========================
## 6B) GENOTYPE-ADJUSTED MATRICES (keep ENV effect)
## =========================
adjust_for_genotype_keep_env <- function(M, geno_factor, env_factor) {
  t(limma::removeBatchEffect(
    t(M),
    batch  = geno_factor,
    design = stats::model.matrix(~ env_factor)
  ))
}

rna_env_adj  <- scale(adjust_for_genotype_keep_env(rna,  geno, env))
prot_env_adj <- scale(adjust_for_genotype_keep_env(prot, geno, env))

## =========================
## 7) DIABLO MODEL (PLAN B)
## =========================
X <- list(
  transcriptome = rna,
  proteome      = prot
)

design <- matrix(c(0, 0.5,
                   0.5, 0),
                 nrow=2, byrow=TRUE)
colnames(design) <- rownames(design) <- names(X)

keepX_best <- list(
  transcriptome = c(50, 50),
  proteome      = c(30, 30)
)

set.seed(1)
diablo <- block.splsda(
  X, Y,
  ncomp  = 2,
  design = design,
  keepX  = keepX_best
)

saveRDS(diablo, file.path(outdir, "DIABLO_model_planB.rds"))

## =========================
## 8) DIABLO FIGURES
## =========================
save_plot_pdf("DIABLO_individuals_planB.pdf", {
  p <- plotIndiv(diablo, ind.names=FALSE, legend=TRUE, title="DIABLO - individuals")
  if (inherits(p,"ggplot")) print(p)
}, width=10, height=8)

save_plot_tiff("DIABLO_individuals_planB.tiff", {
  p <- plotIndiv(diablo, ind.names=FALSE, legend=TRUE, title="DIABLO - individuals")
  if (inherits(p,"ggplot")) print(p)
}, width_in=10, height_in=8)

save_plot_pdf("DIABLO_variables_planB.pdf", {
  p <- plotVar(diablo, var.names=FALSE, title="DIABLO - variables")
  if (inherits(p,"ggplot")) print(p)
}, width=10, height=8)

save_plot_tiff("DIABLO_variables_planB.tiff", {
  p <- plotVar(diablo, var.names=FALSE, title="DIABLO - variables")
  if (inherits(p,"ggplot")) print(p)
}, width_in=10, height_in=8)

save_plot_pdf("DIABLO_block_correlations_planB.pdf", {
  p <- plotDiablo(diablo, ncomp=2)
  if (inherits(p,"ggplot")) print(p)
}, width=12, height=6.5)

save_plot_tiff("DIABLO_block_correlations_planB.tiff", {
  p <- plotDiablo(diablo, ncomp=2)
  if (inherits(p,"ggplot")) print(p)
}, width_in=12, height_in=6.5)

save_plot_pdf("DIABLO_arrow_planB.pdf", {
  p <- plotArrow(diablo, ind.names=FALSE, legend=TRUE, title="DIABLO - arrow (comp1 -> comp2)")
  if (inherits(p,"ggplot")) print(p)
}, width=10, height=8)

save_plot_tiff("DIABLO_arrow_planB.tiff", {
  p <- plotArrow(diablo, ind.names=FALSE, legend=TRUE, title="DIABLO - arrow (comp1 -> comp2)")
  if (inherits(p,"ggplot")) print(p)
}, width_in=10, height_in=8)

save_plot_pdf("DIABLO_cim_comp1_planB.pdf", { cimDiablo(diablo, comp=1, margin=c(10,10)) },
              width=14, height=9)
save_plot_tiff("DIABLO_cim_comp1_planB.tiff", { cimDiablo(diablo, comp=1, margin=c(10,10)) },
               width_in=14, height_in=9, dpi=600)

save_plot_pdf("DIABLO_cim_comp2_planB.pdf", { cimDiablo(diablo, comp=2, margin=c(10,10)) },
              width=14, height=9)
save_plot_tiff("DIABLO_cim_comp2_planB.tiff", { cimDiablo(diablo, comp=2, margin=c(10,10)) },
               width_in=14, height_in=9, dpi=600)

cutoff_circos <- 0.7
save_plot_pdf("DIABLO_circos_comp1_cutoff0p7.pdf", {
  circosPlot(diablo, comp=1, cutoff=cutoff_circos, line=TRUE,
             color.blocks=c("dodgerblue3","firebrick2"),
             color.cor=c("chocolate3","grey20"),
             size.labels=1.5)
  title(main=paste0("DIABLO circos - Comp1 (cutoff ", cutoff_circos, ")"), cex.main=1.2)
}, width=14, height=10)

save_plot_tiff("DIABLO_circos_comp1_cutoff0p7.tiff", {
  circosPlot(diablo, comp=1, cutoff=cutoff_circos, line=TRUE,
             color.blocks=c("dodgerblue3","firebrick2"),
             color.cor=c("chocolate3","grey20"),
             size.labels=1.5)
  title(main=paste0("DIABLO circos - Comp1 (cutoff ", cutoff_circos, ")"), cex.main=1.2)
}, width_in=14, height_in=10, dpi=600)

## =========================
## 9) EXPORT TABLES (CSV + Excel)
## =========================
sel_all <- rbind(
  get_selected_features(diablo, "transcriptome", 1),
  get_selected_features(diablo, "proteome",      1),
  get_selected_features(diablo, "transcriptome", 2),
  get_selected_features(diablo, "proteome",      2)
)

if (nrow(sel_all) > 0) {
  sel_all$abs_loading <- abs(sel_all$loading)
  sel_all <- sel_all[order(sel_all$block, sel_all$comp, -sel_all$abs_loading), ]
  sel_all$abs_loading <- NULL
}
write_csv(sel_all, "DIABLO_selected_features_loadings_all.csv")

cc1 <- make_cross_cor(diablo, X, comp=1)
cc2 <- make_cross_cor(diablo, X, comp=2)
if (nrow(cc1$long) > 0) write_csv(cc1$long, "DIABLO_crosscor_long_comp1.csv")
if (nrow(cc2$long) > 0) write_csv(cc2$long, "DIABLO_crosscor_long_comp2.csv")

xlsx_path <- file.path(outdir, "DIABLO_outputs_all_tables.xlsx")
wb <- createWorkbook()
add_sheet <- function(name, df) { addWorksheet(wb, name); writeData(wb, sheet=name, x=df) }
add_sheet("selected_features_all", sel_all)
if (nrow(cc1$long) > 0) add_sheet("crosscor_long_comp1", cc1$long)
if (nrow(cc2$long) > 0) add_sheet("crosscor_long_comp2", cc2$long)
saveWorkbook(wb, xlsx_path, overwrite=TRUE)

## =========================
## 10) HEATMAPS — selected DIABLO features (no sample clustering)
## =========================
sel_rna_c1  <- unique(get_selected_features(diablo, "transcriptome", 1)$feature); sel_rna_c1  <- sel_rna_c1[!is.na(sel_rna_c1)]
sel_prot_c1 <- unique(get_selected_features(diablo, "proteome",      1)$feature); sel_prot_c1 <- sel_prot_c1[!is.na(sel_prot_c1)]
sel_rna_c2  <- unique(get_selected_features(diablo, "transcriptome", 2)$feature); sel_rna_c2  <- sel_rna_c2[!is.na(sel_rna_c2)]
sel_prot_c2 <- unique(get_selected_features(diablo, "proteome",      2)$feature); sel_prot_c2 <- sel_prot_c2[!is.na(sel_prot_c2)]

rna_c1_mat  <- X$transcriptome[, intersect(sel_rna_c1,  colnames(X$transcriptome)), drop=FALSE]
prot_c1_mat <- X$proteome[,      intersect(sel_prot_c1, colnames(X$proteome)),      drop=FALSE]
rna_c2_mat  <- X$transcriptome[, intersect(sel_rna_c2,  colnames(X$transcriptome)), drop=FALSE]
prot_c2_mat <- X$proteome[,      intersect(sel_prot_c2, colnames(X$proteome)),      drop=FALSE]

make_heatmap_nocluster(rna_c1_mat,  group=Y, filename_prefix="HEATMAP_RNA_selected_comp1_noSampleCluster",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(prot_c1_mat, group=Y, filename_prefix="HEATMAP_PROT_selected_comp1_noSampleCluster",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(rna_c2_mat,  group=Y, filename_prefix="HEATMAP_RNA_selected_comp2_noSampleCluster",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(prot_c2_mat, group=Y, filename_prefix="HEATMAP_PROT_selected_comp2_noSampleCluster",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)

## ENV-focused heatmaps (genotype removed)
rna_env_c1_mat  <- rna_env_adj[,  intersect(sel_rna_c1,  colnames(rna_env_adj)),  drop=FALSE]
prot_env_c1_mat <- prot_env_adj[, intersect(sel_prot_c1, colnames(prot_env_adj)), drop=FALSE]
rna_env_c2_mat  <- rna_env_adj[,  intersect(sel_rna_c2,  colnames(rna_env_adj)),  drop=FALSE]
prot_env_c2_mat <- prot_env_adj[, intersect(sel_prot_c2, colnames(prot_env_adj)), drop=FALSE]

make_heatmap_nocluster(rna_env_c1_mat,  group=Y_env, filename_prefix="HEATMAP_ENVFOCUSED_RNA_selected_comp1_genoRemoved",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(prot_env_c1_mat, group=Y_env, filename_prefix="HEATMAP_ENVFOCUSED_PROT_selected_comp1_genoRemoved",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(rna_env_c2_mat,  group=Y_env, filename_prefix="HEATMAP_ENVFOCUSED_RNA_selected_comp2_genoRemoved",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
make_heatmap_nocluster(prot_env_c2_mat, group=Y_env, filename_prefix="HEATMAP_ENVFOCUSED_PROT_selected_comp2_genoRemoved",
                       scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)

## =========================
## 11) CIRBP AXIS — Fig6B / Fig6C / ExtFig6-3
## =========================
find_feature <- function(cols, target="CIRBP") {
  if (target %in% cols) return(target)
  hit <- cols[grepl(paste0("^", target, "$"), cols, ignore.case=TRUE)]
  if (length(hit)) return(hit[1])
  hit <- cols[grepl(target, cols, ignore.case=TRUE)]
  if (length(hit)) return(hit[1])
  NA_character_
}

spearman_cor <- function(x, y) suppressWarnings(cor(x, y, method="spearman", use="pairwise.complete.obs"))
zscore_cols <- function(M) scale(M)

make_cirbp_module <- function(M, cirbp_col, top_n=40, min_abs_r=0.60) {
  r <- apply(M, 2, function(v) spearman_cor(v, M[, cirbp_col]))
  r <- r[names(r) != cirbp_col]
  r <- r[order(abs(r), decreasing=TRUE)]
  r <- r[abs(r) >= min_abs_r]
  feats <- head(names(r), top_n)
  list(features=feats, r=r[feats])
}

module_score <- function(M, features, weights=NULL) {
  features <- intersect(features, colnames(M))
  if (length(features) == 0) return(rep(NA_real_, nrow(M)))
  Z <- zscore_cols(M[, features, drop=FALSE])
  if (is.null(weights)) weights <- rep(1, length(features))
  weights <- weights[seq_along(features)]
  weights <- weights / sum(abs(weights))
  as.numeric(Z %*% weights)
}

## locate CIRBP
cirbp_rna  <- find_feature(colnames(X$transcriptome), "CIRBP")
cirbp_prot <- find_feature(colnames(X$proteome),      "CIRBP")
cat("\n[CIRBP lookup]\nRNA:  ", cirbp_rna, "\nPROT: ", cirbp_prot, "\n\n", sep="")

rna_mod  <- if (!is.na(cirbp_rna))  make_cirbp_module(X$transcriptome, cirbp_rna,  top_n=50, min_abs_r=0.60) else NULL
prot_mod <- if (!is.na(cirbp_prot)) make_cirbp_module(X$proteome,      cirbp_prot, top_n=30, min_abs_r=0.60) else NULL

scores <- data.frame(sample=meta$sample, stringsAsFactors=FALSE)
rownames(scores) <- scores$sample

if (!is.null(rna_mod)) {
  w <- sign(rna_mod$r); w[is.na(w)] <- 1
  scores$CIRBP_RNA_module <- module_score(X$transcriptome, rna_mod$features, weights=w)
  scores$CIRBP_RNA_level  <- X$transcriptome[, cirbp_rna]
}
if (!is.null(prot_mod)) {
  w <- sign(prot_mod$r); w[is.na(w)] <- 1
  scores$CIRBP_PROT_module <- module_score(X$proteome, prot_mod$features, weights=w)
  scores$CIRBP_PROT_level  <- X$proteome[, cirbp_prot]
}

avail <- intersect(c("CIRBP_RNA_module","CIRBP_PROT_module"), colnames(scores))
scores$CIRBP_multiomic_module <- rowMeans(scores[, avail, drop=FALSE], na.rm=TRUE)

scores <- cbind(scores, meta[rownames(scores), c("genotype","environment","group")])
scores$genotype    <- factor(scores$genotype)
scores$environment <- factor(scores$environment, levels=c("GROUND","ISS"))

## Fig 6B
p6B <- ggplot(scores, aes(x=environment, y=CIRBP_multiomic_module)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.12, size=1.8) +
  facet_wrap(~ genotype, nrow=1) +
  labs(title="Fig 6B | CIRBP-linked multi-omic module score", x=NULL, y="Module score (z)") +
  theme_bw()

save_plot_pdf("Fig6B_CIRBP_module_byENV_facetGenotype.pdf", { print(p6B) }, width=12, height=5)
save_plot_tiff("Fig6B_CIRBP_module_byENV_facetGenotype.tiff", { print(p6B) }, width_in=12, height_in=5, dpi=600)
write_csv(scores, "Fig6B_CIRBP_module_scores_allSamples.csv")

## Ext Fig 6-3: CIRBP module heatmaps
if (!is.null(rna_mod) && length(rna_mod$features) > 0) {
  mat_rna_mod <- X$transcriptome[, intersect(rna_mod$features, colnames(X$transcriptome)), drop=FALSE]
  make_heatmap_nocluster(mat_rna_mod, group=Y,
                         filename_prefix="ExtFig6_3A_CIRBP_RNA_module_heatmap",
                         scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
  write_csv(data.frame(feature=colnames(mat_rna_mod), block="transcriptome"), "ExtFig6_3A_CIRBP_RNA_module_features.csv")
}
if (!is.null(prot_mod) && length(prot_mod$features) > 0) {
  mat_prot_mod <- X$proteome[, intersect(prot_mod$features, colnames(X$proteome)), drop=FALSE]
  make_heatmap_nocluster(mat_prot_mod, group=Y,
                         filename_prefix="ExtFig6_3B_CIRBP_PROT_module_heatmap",
                         scale_rows=TRUE, show_rownames=FALSE, cluster_features=TRUE)
  write_csv(data.frame(feature=colnames(mat_prot_mod), block="proteome"), "ExtFig6_3B_CIRBP_PROT_module_features.csv")
}

## =========================
## 12) PATHWAYS — ONLY THE 5 MANUAL SETS YOU CHOSE
## =========================
pathways_list <- list(
  
  "Cellular Senescence" = c(
    "TFDP1","CDKN1B","TNRC6B","H2BC12L","MAPK7","UBA52","TNRC6C","HMGA1","LMNB1","ASF1A",
    "TXN","AGO4","ANAPC16","TNRC6A","KDM6B","CDKN2D","RELA","H2AC6","MOV10","H3C1",
    "CBX4","MINK1","UBE2C","RPS27A","PHC3","H1-0","MAPKAPK2","H2BC26","H4C4","H2AX",
    "TNIK","RBBP4","MAPK14","JUN","ANAPC7","RBBP7","STAT3","TINF2","HMGA2","UBE2S",
    "H2BC15","KAT5","CBX8","ACD","UBE2E1","RPS6KA3","MAPK3","MAPKAPK5","MAPK1","RING1",
    "CBX6","H1-3","H1-5","RNF2","EHMT1","MAPK11","MAP4K4","MAP2K7","MAPK10","TERF2",
    "MAPK9","UBE2D1","MAP2K4","AGO3","AGO1","FZR1","RAD50","SUZ12","TERF2IP","CDC23",
    "EZH2","RPS6KA1","PHC2","MAPK8","CDC16","CDK4","PHC1","CABIN1","TERF1","ANAPC5",
    "MRE11","CDK6","ANAPC2","ATM","MAPKAPK3","MAP2K6","SCMH1","EHMT2","ANAPC1","H1-1",
    "RB1","EED","RPS6KA2","ANAPC4","CDC27","HIRA","NBN","CDK2","SP1","POT1",
    "CBX2","TP53","MAP2K3","NFKB1","ANAPC10","CEBPB","IGFBP7","E2F3","TFDP2"
  ),
  
  "NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION" = c(
    "GRIN1","GRIN2A","GRIN2B","GRIN2C","GRIN2D",
    "GRIA1","GRIA2","GRIA3","GRIA4",
    "GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6",
    "GABRB1","GABRB2","GABRB3","GABRG2","GABRG3",
    "DLG1","DLG2","DLG3","DLG4",
    "CAMK2A","CAMK2B","CAMK2D","CAMK2G",
    "CREB1","PRKCA","PRKCB","PRKCG","MAPK1","MAPK3"
  ),
  
  "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES" = c(
    "PPP1R13B","TP53AIP1","ZNF420","CASP2","BCL2L14","PIDD1","CREBBP","APAF1","TP63",
    "RABGGTA","BCL6","NLRC4","CASP10","CASP1","PRELID3A","BNIP3L","STEAP3","IGFBP3",
    "PRELID1","TRIAP1","CRADD","BID","TP73","TNFRSF10A","ATM","TNFRSF10D","RABGGTB",
    "PERP","TMEM219","BIRC5","BBC3","TNFRSF10C","NDRG1","CASP6","TP53I3","CHM",
    "AIFM2","TP53INP1","BAX","TP53","TNFRSF10B","TP53BP2","PMAIP1","FAS"
  ),
  
  "METABOLISM_OF_RNA" = c(
    "LAS1L","RBM5","UPF1","POLR2J","GCFC2","ELAC2","TSR3","PSMB1","RPS20","UTP18",
    "PTBP1","PSMC4","ZC3H3","STEEP1","MNAT1","AQR","TFB1M","NUP160","ZCCHC8","RPL26L1",
    "NSUN2","METTL1","EDC4","MTREX","PSMA4","TRIT1","GEMIN8","YTHDC2","POLR2B","TPR",
    "THOC3","NOP58","ZC3H11A","RIOK2","NDC1","CCAR1","SNRNP40","EIF4B","RPL18","U2AF2",
    "DDX20","WDR3","WDR18","ADAT1","YBX1","ELAVL1","TRMT11","THUMPD1","DHX8","GPKOW",
    "FTSJ1","NUP133","SMG6","PABPC1","RPL31","BUD23","FAM50A","CPSF1","HSD17B10","TRNT1"
  ),
  
  "OXIDATIVE_STRESS_RESPONSE" = c(
    "NFE2L2","KEAP1","SOD1","SOD2","CAT","GPX1","GPX4","TXN","TXNRD1",
    "PRDX1","PRDX2","PRDX3","PRDX5","PRDX6","HMOX1","HMOX2","NQO1",
    "GSR","GCLC","GCLM"
  ),
  
  "AUTOPHAGY_LYSOSOME" = c(
    "ATG3","ATG5","ATG7","ATG9A","ATG12","ATG13","BECN1",
    "ULK1","ULK2","MAP1LC3A","MAP1LC3B","SQSTM1",
    "LAMP1","LAMP2","CTSB","CTSD","CTSL"
  ),
  
  "SYNAPTIC_VESICLE_CYCLE" = c(
    "SYN1","SYN2","SYN3",
    "SNAP25","STX1A","STX1B",
    "VAMP1","VAMP2",
    "SLC17A7","SLC17A6",
    "SYT1","SYT2",
    "DNM1","DNM2",
    "NSF"
  )
  
)

compute_proxy_score <- function(M, genes, min_genes=6) {
  genes_present <- intersect(unique(genes), colnames(M))
  if (length(genes_present) < min_genes) {
    return(list(score=rep(NA_real_, nrow(M)),
                n=length(genes_present),
                genes_present=genes_present))
  }
  Z <- scale(M[, genes_present, drop=FALSE])
  list(score=rowMeans(Z, na.rm=TRUE),
       n=length(genes_present),
       genes_present=genes_present)
}

run_pathway_pipeline <- function(pathway_name, genes, min_genes=6) {
  
  tag  <- sanitize_tag(pathway_name)
  ycol <- paste0("PW_", tag)
  
  M <- X$transcriptome[rownames(scores), , drop=FALSE]  # align to scores order
  res <- compute_proxy_score(M, genes, min_genes=min_genes)
  
  if (all(is.na(res$score))) {
    warning("Skipping ", pathway_name, " (only ", res$n, " genes present).")
    return(NULL)
  }
  
  df <- scores
  df[[ycol]] <- res$score
  
  ## save genes used
  write_csv(data.frame(pathway=pathway_name, gene=res$genes_present, stringsAsFactors=FALSE),
            paste0("GENES_USED_", tag, ".csv"))
  
  df2 <- df %>%
    filter(!is.na(.data[[ycol]]),
           !is.na(genotype), !is.na(environment),
           !is.na(CIRBP_multiomic_module))
  
  ## ANOVA geno × env
  model <- lm(as.formula(paste0(ycol, " ~ genotype * environment")), data=df2)
  aov_tab <- anova(model)
  
  terms <- setdiff(rownames(aov_tab), "Residuals")
  ss_res <- as.numeric(aov_tab["Residuals","Sum Sq"])
  ss_terms <- as.numeric(aov_tab[terms,"Sum Sq"])
  names(ss_terms) <- terms
  partial_eta2 <- ss_terms / (ss_terms + ss_res)
  
  aov_out <- data.frame(
    pathway=pathway_name,
    genes_in_set=length(unique(genes)),
    genes_present=res$n,
    term=terms,
    p_value=as.numeric(aov_tab[terms,"Pr(>F)"]),
    partial_eta2=as.numeric(partial_eta2),
    stringsAsFactors=FALSE
  )
  write_csv(aov_out, paste0("ANOVA_", tag, "_genotype_environment_effectSizes.csv"))
  
  ## posthoc env overall + within genotype
  emm_env <- emmeans(model, ~ environment)
  post_overall <- as.data.frame(pairs(emm_env, adjust="BH"))
  post_overall$pathway <- pathway_name
  write_csv(post_overall, paste0("POSTHOC_", tag, "_ISS_vs_GROUND_overall.csv"))
  
  emm_env_by_gt <- emmeans(model, ~ environment | genotype)
  post_by_gt <- as.data.frame(pairs(emm_env_by_gt, adjust="BH"))
  post_by_gt$pathway <- pathway_name
  write_csv(post_by_gt, paste0("POSTHOC_", tag, "_ISS_vs_GROUND_withinGenotype.csv"))
  
  ## correlation with CIRBP module
  cor_global <- spearman_cor(df2$CIRBP_multiomic_module, df2[[ycol]])
  cor_ground <- spearman_cor(df2$CIRBP_multiomic_module[df2$environment=="GROUND"],
                             df2[[ycol]][df2$environment=="GROUND"])
  cor_iss <- spearman_cor(df2$CIRBP_multiomic_module[df2$environment=="ISS"],
                          df2[[ycol]][df2$environment=="ISS"])
  
  cor_out <- data.frame(
    pathway=pathway_name,
    genes_present=res$n,
    spearman_r_global=cor_global,
    spearman_r_ground=cor_ground,
    spearman_r_iss=cor_iss,
    stringsAsFactors=FALSE
  )
  write_csv(cor_out, paste0("CORR_CIRBP_vs_", tag, ".csv"))
  
  ## plots: boxplot by env facet genotype
  p1 <- ggplot(df2, aes(x=environment, y=.data[[ycol]], color=environment)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.12, size=1.8) +
    facet_wrap(~ genotype, nrow=1) +
    theme_bw() +
    labs(title=paste0(pathway_name, " | ISS vs GROUND (proxy z)"),
         x=NULL, y=paste0(pathway_name, " activity (proxy z)")) +
    theme(legend.position="none")
  
  save_plot_pdf(paste0("Fig_", tag, "_ISS_vs_GROUND_facetGenotype.pdf"), { print(p1) }, width=12, height=5)
  save_plot_tiff(paste0("Fig_", tag, "_ISS_vs_GROUND_facetGenotype.tiff"), { print(p1) }, width_in=12, height_in=5, dpi=600)
  
  ## interaction plot (estimated means)
  emm_int <- emmeans(model, ~ environment | genotype)
  emm_int_df <- as.data.frame(emm_int)
  
  p2 <- ggplot(emm_int_df, aes(x=environment, y=emmean, group=genotype)) +
    geom_line() +
    geom_point(size=2.5) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.08) +
    theme_bw() +
    labs(title=paste0("Interaction | ", pathway_name, " (estimated means)"),
         x="Environment", y="Estimated mean (proxy z)")
  
  save_plot_pdf(paste0("Fig_", tag, "_interactionPlot_emmeans.pdf"), { print(p2) }, width=7.5, height=5.5)
  save_plot_tiff(paste0("Fig_", tag, "_interactionPlot_emmeans.tiff"), { print(p2) }, width_in=7.5, height_in=5.5, dpi=600)
  
  ## scatter vs CIRBP (Fig 6C-style)
  p3 <- ggplot(df2, aes(x=.data[[ycol]], y=CIRBP_multiomic_module, color=environment, shape=genotype)) +
    geom_point(size=2) +
    theme_bw() +
    labs(title=paste0("Fig 6C | CIRBP module vs ", pathway_name),
         x=paste0(pathway_name, " activity (proxy z)"),
         y="CIRBP multi-omic module (z)")
  
  save_plot_pdf(paste0("Fig6C_CIRBP_vs_", tag, "_scatter.pdf"), { print(p3) }, width=7.5, height=5.5)
  save_plot_tiff(paste0("Fig6C_CIRBP_vs_", tag, "_scatter.tiff"), { print(p3) }, width_in=7.5, height_in=5.5, dpi=600)
  
  cat("\n✅ Done pathway: ", pathway_name,
      " | genes present: ", res$n,
      " | r(CIRBP): ", round(cor_global, 3),
      "\n", sep="")
  
  invisible(list(aov=aov_out, corr=cor_out))
}

## Run the 5 pathways
  res_list <- lapply(names(pathways_list), function(nm) {
    
    min_g <- dplyr::case_when(
      nm == "Cellular Senescence" ~ 10,
      nm == "METABOLISM_OF_RNA" ~ 15,
      nm == "NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION" ~ 6,
      nm %in% c(
        "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
        "OXIDATIVE_STRESS_RESPONSE",
        "AUTOPHAGY_LYSOSOME",
        "SYNAPTIC_VESICLE_CYCLE"
      ) ~ 4,
      TRUE ~ 6
    )
    
    run_pathway_pipeline(
      pathway_name = nm,
      genes = pathways_list[[nm]],
      min_genes = min_g
    )
    
  })
  
  names(res_list) <- names(pathways_list)
names(res_list) <- names(pathways_list)

## Summary correlation table
summ <- do.call(rbind, lapply(res_list, function(x) if (is.null(x)) NULL else x$corr))
if (!is.null(summ) && nrow(summ) > 0) {
  write_csv(summ, "SUMMARY_manualPathways_CIRBP_correlations.csv")
  print(summ)
}

## =========================
## 14) CELL-TYPE MARKER ANALYSIS (RNA + PROTEOME)
## based on marker table provided by user
## =========================

celltype_markers <- list(
  "radial_glia" = c("SFRP1","ZFP36L1","ID4","HMGB2","HES1","TTYH1","SOX3","SYNE2","LHX2","HMGA2",
                    "SOX2","C1orf61","B3GAT2","GLI3","CENPF","TOP2A","UBE2C","CKB","LINC01551",
                    "HIST1H3D","HIST1H4C","TCF7L1","PCLAF","AMBN","GAS1","CREB5","HIST1H1B","HES5",
                    "EMX2","SOX9","GNAI2","IFI44L","NAP1L1","PTTG1","QKI","CKS2","DOK5","PAX6","HEY1",
                    "NUSAP1","TMEM161B-AS1","HIST1H3B","SMC4","DMRTA2","NASP","LINC00461","ASPM",
                    "HIST1H1C","CRYGD","H2AFX"),
  
  "glutamatergic_neuron_intermediate_progenitor" = c("NEUROD6","NEUROD2","SSTR2","TBR1","ZBTB18","NHLH1",
                                                     "IGFBPL1","NRN1","RTN1","THSD7A","NRXN1","BHLHE22",
                                                     "CALB2","KHDRBS3","CCSAP","PDE1A","NEUROD1","NPTX1",
                                                     "NXPH4","NTS","NEUROG2","OLFM1","PRDM8","CORO2B",
                                                     "TP53I11","ZFPM2","PCDH9","NELL2","SRRM4","SCG3",
                                                     "DCC","EPB41L3","SLC17A7","ST18","EOMES","NSG2",
                                                     "EMX1","CAP2","SYT4","SPINK5","NSMF","ANK3","MYT1L",
                                                     "FSTL5","EBF2","CELF4","B3GAT1","EPHA5","NHLH2","DLL3"),
  
  "Cajal_Retzius_cell" = c("RELN","LHX1","NHLH2","PCP4","RSPO3","NR2F2","PGF","MAB21L1","CDKN1A",
                           "LHX5_AS1","SEMA6A","EBF3","LHX1_DT","TBR1","PIDD1","EMX2","MYT1L","OLFM1",
                           "THSD7A","RTN1","INA","SPINK5","NHLH1","CXCR4","TAGLN2","LHX5","TP73",
                           "WNT7B","NRN1","RGS16","IGSF9","DCC","FNDC5","SVEP1","SIPA1L2","ELAVL2",
                           "LMO2","SCG3","SAMD3","CACNA2D1","CDK5R1","KLF7","DMRTA2","CPNE4","SNAP25",
                           "GNG3","CADM2","BHLHE22","CABP7"),
  
  "GABAergic_progenitor_early" = c("HIST1H4C","TOP2A","HIST1H3B","HIST1H3D","HIST1H2AG","NUSAP1",
                                   "HIST1H1B","HIST2H2AC","GAD2","CENPF","MKI67","HIST1H1C","SMC4",
                                   "HMGN2","HMGB2","TUBA1B","DLX2","HIST1H4D","DLX5","ASCL1","DLX6_AS1",
                                   "PBK","ASPM","PRC1","ZEB2","CDK1","TPX2","DLX1","MAD2L1","NPY",
                                   "HMGB1","CKS1B","HIST1H4E","BIRC5","PFN2","HIST1H1E","SHTN1","DLGAP5",
                                   "AURKB","UBE2C","KIF11","C11orf96","UBE2S","GTSE1","CENPE","CKS2",
                                   "FBXO5","CKAP2","H2AFX","KNL1"),
  
  "GABAergic_progenitor_late" = c("DLX2","PANTR1","DLX1","CENPF","TOP2A","DLX5","ASCL1","ZEB2","GAD2",
                                  "SLC1A3","FAM181B","NUSAP1","HIST1H1B","DLX6_AS1","MKI67","BIRC5","NPY",
                                  "GSX2","UBE2C","ASPM","LIX1","PBK","SCRG1","MT3","RORB","CDO1","ADGRV1",
                                  "LINC01896","HOPX","PRC1","ETV1","TPX2","RAB31","TMEM158","TNC","CDK1",
                                  "SCGN","SMOC1","CENPE","ATAD2","NUF2","NES","CDC20","STK39","MOXD1",
                                  "PBX3","PCLAF","CDKN3","NTM"),
  
  "GABAergic_neuron" = c("SIX3","DLX6_AS1","SST","PBX3","TAC1","RUNX1T1","STMN2","NSG2","RTN1","DLX5",
                         "GRIA2","INA","MAPT","TTC9B","AC004158.1","DLX6","GAD2","NRGN","SCG5","MYT1L",
                         "SCG3","RBFOX1","SLC8A1","RAB3A","CDK5R1","GRIN2B","CELF4","FOXP2","GABRB3",
                         "DCLK1","MTSS1","NRXN3","MEF2C","TMEM35A","SERTAD4","SP9","ISLR2","DPF1",
                         "KIF21B","ATP1A3","SCN3B","SLITRK6","DNM3","CXXC4","PCDH10","ACTL6B","SCN3A",
                         "LINGO1","CELF3"),
  
  "caudal_neuron" = c("CRABP1","NEFL","NEFM","ROBO3","CNTN2","HOXA5","RGMB","HOXB5","LINC02381","LHX1",
                      "POU4F1","TFAP2A","HOXB3","LINC00682","ZIC1","PHOX2B","HOXB4","HOTAIRM1","CACNA2D1",
                      "MAB21L2","HOXB2","UNCX","NRN1","HOXA3","TFAP2B","EBF3","GDAP1L1","NHLH1","ELAVL2",
                      "DCC","HOXA4","EBF1","KLHL35","HOXC4","AC008522.1","LAMP5","NHLH2","PRPH","FNDC5",
                      "POU3F1","IRX2","POU3F2","UBE2QL1","SCRT2","SCG3","NRCAM","SNAP25","CBLN1","MYT1"),
  
  "choroid_plexus" = c("TTR","TPBG","TRPM3","FAM89A","RSPO3","RSPO2","WLS","IGFBP7","SOX2_OT","ID1",
                       "ID3","NEAT1","HTR2C","ZFP36L2","SULF1","PLS3","OTX2","TRIM67","PCP4","CCK","EMX2",
                       "CA2","PRTG","CLDN5","CLU","PMCH","CRABP2","MAF","BMP7","MSX1","FOS","CACHD1",
                       "EGR1","RMST","CXCL14","SLC7A8","ATP1A2","ANKRD37","ZIC2","CDO1","SPRY1","PIFO",
                       "NR2F2","LGI1","TSHZ2","DCT","EMX2OS","KLHDC8B","GJA1","SDC2"),
  
  "CNS_fibroblast_1" = c("COL3A1","COL1A1","DCN","SPARC","COL1A2","LUM","MGP","DLK1","FRZB","POSTN",
                         "H19","S100A10","APOE","SLC9A3R1","IFITM3","IGFBP3","KRT19","FN1","KRT8","COL5A1",
                         "TIMP1","PLTP","KRT18","COL6A2","ANXA2","OLFML3","COL6A1","TMEM88","IL6ST","CXCL14",
                         "S100A11","COL15A1","PDPN","B2M","HAND1","PODXL","COL5A2","SERPINH1","AQP1","RRBP1",
                         "EGFL6","SLC40A1","CHN1","TFPI2","EMP3","COL21A1","PITX1","BST2","CDKN1C"),
  
  "CNS_fibroblast_2" = c("GABRP","S100A11","KRT19","KRT18","KRT8","FN1","IGFBP3","S100A10","ANXA1","CRABP2",
                         "HPGD","NPY","ANXA2","CLDN6","MT1E","EPCAM","IGFBP7","EMP2","ANXA3","TGFBI","FOS",
                         "AHNAK","EPAS1","ANKRD1","MT2A","CCN1","LGALS3","CTHRC1","COL12A1","APOE","SLC2A3",
                         "ARHGAP29","PRSS23","LY6E","CCDC80","SCUBE3","L1TD1","KCNMA1","TNC","TSTD1","TAGLN",
                         "FSTL1","CXCL12","TPM2","VTCN1","SPATS2L","IGDCC3","CLDN4","PERP","S100A13"),
  
  "neural_crest_early" = c("HES5","ZIC1","MSX1","NES","ATP1A2","MIR99AHG","CCND1","NTRK2","ID1","ID3",
                           "RFX4","CRNDE","OLIG3","AC092957.1","GRID2","LINC02381","NCKAP5","AL139246.5",
                           "TOP2A","CDC25B","BOC","IRX2","MASP1","SLC1A3","WLS","CXCR4","ARL4A","FZD10",
                           "ZIC4","RND2","TPBG","GYPC","PAX3","MEIS1","IRX1","DLL1","ZNF503","EDNRB",
                           "HOXB3","MKI67","PCDH18","LRIG1","MGST1","EPHA3","NRARP","JAG1","NAV2","HOXB2",
                           "AL359091.1","APCDD1"),
  
  "neural_crest_middle" = c("LGALS1","TWIST1","PRRX1","SHOX2","GPC3","DNM3OS","PDGFRA","ID1","CRABP2","ID3",
                            "KCTD12","S100A11","NR2F2","PCDH18","ID2","CDC42EP5","TPM2","LRRC17","LIMA1",
                            "MECOM","EPB41L2","COL1A2","TES","CCDC80","HOXA10","TGFB2","MIR100HG","CXCL12",
                            "HMGA2","MIR99AHG","CDH11","IGF2","NRP2","SFRP2","MFAP2","SIX1","IL11RA",
                            "ZFHX4_AS1","PCOLCE","C1orf54","HMCN1","PCDH10","TBX3","ITM2C","FGFR1","RAB34",
                            "IRX3","EVA1B","ARHGAP28","CAST"),
  
  "neural_crest_late" = c("S100A6","S100B","OLFML2A","AL139393.3","PLP1","S100A10","FN1","BCHE","S100A4","SOX10",
                          "ADAM12","S100A11","IGFBP7","HIST1H1A","LGALS1","MMP2","ERBB3","TFAP2A","COL4A2",
                          "LAMC1","ARPC1B","CDH19","SLITRK6","VIM","HIST1H3D","CRYL1","TAC1","RXRG","BCAN",
                          "PRSS23","DST","ITM2C","ZEB2","KCTD12","COL4A1","HMG20B","LMNA","CLIC1","NR2F2",
                          "SLC25A5","AHR","RHOC","EDNRB","MMP17","PMP22","MIA","FAM210B","CST3","APOD","GAS7"),
  
  "endothelial_cell" = c("H19","KDR","TFPI","CARTPT","CAV1","CD93","PECAM1","GJA4","GNG11","CDH5","CLDN5",
                         "ICAM2","IGFBP4","PRND","FLT1","ESAM","COL4A1","IFITM3","RNASE1","TIMP3","PRSS23",
                         "CD34","PLVAP","RAMP2","CRIM1","FN1","EGFL7","ACSM3","TIE1","VAMP5","CALCRL","S100A16",
                         "CAVIN1","IGF2","HLA_E","HSPG2","IFITM2","MECOM","ADGRL4","COL4A2","NRP1","BST2",
                         "LIMS1","MYH9","PLXND1","IFI16","ARHGAP29","ECSCR","EPAS1","CLEC14A")
)

## helpers to harmonize names a bit
normalize_features <- function(x) {
  x <- gsub("-", "_", x)
  x <- gsub("\\.", "_", x)
  toupper(x)
}

match_features_caseflex <- function(targets, cols) {
  targets2 <- normalize_features(targets)
  cols2    <- normalize_features(cols)
  m <- match(targets2, cols2)
  cols[unique(stats::na.omit(m))]
}

compute_celltype_score <- function(M, genes, min_genes=5) {
  genes_present <- match_features_caseflex(genes, colnames(M))
  if (length(genes_present) < min_genes) {
    return(list(score=rep(NA_real_, nrow(M)),
                n=length(genes_present),
                genes_present=genes_present))
  }
  Z <- scale(M[, genes_present, drop=FALSE])
  list(score=rowMeans(Z, na.rm=TRUE),
       n=length(genes_present),
       genes_present=genes_present)
}

## output containers
celltype_scores <- data.frame(sample=meta$sample, stringsAsFactors=FALSE)
rownames(celltype_scores) <- celltype_scores$sample
celltype_scores <- cbind(celltype_scores, meta[rownames(celltype_scores), c("genotype","environment","group")])

celltype_summary <- list()
celltype_anova_list <- list()

## loop over cell types
for (ct in names(celltype_markers)) {
  
  tag <- sanitize_tag(ct)
  
  ## RNA score
  res_rna <- compute_celltype_score(X$transcriptome[rownames(celltype_scores), , drop=FALSE],
                                    celltype_markers[[ct]],
                                    min_genes=5)
  
  ## PROT score
  res_prot <- compute_celltype_score(X$proteome[rownames(celltype_scores), , drop=FALSE],
                                     celltype_markers[[ct]],
                                     min_genes=3)
  
  ## save marker presence tables
  write_csv(
    data.frame(cell_type=ct,
               omic="RNA",
               marker_present=res_rna$genes_present,
               stringsAsFactors=FALSE),
    paste0("CELLTYPE_", tag, "_RNA_markers_present.csv")
  )
  
  write_csv(
    data.frame(cell_type=ct,
               omic="PROTEOME",
               marker_present=res_prot$genes_present,
               stringsAsFactors=FALSE),
    paste0("CELLTYPE_", tag, "_PROT_markers_present.csv")
  )
  
  ## add scores
  rna_col  <- paste0("CT_", tag, "_RNA")
  prot_col <- paste0("CT_", tag, "_PROT")
  comb_col <- paste0("CT_", tag, "_COMBINED")
  
  celltype_scores[[rna_col]]  <- res_rna$score
  celltype_scores[[prot_col]] <- res_prot$score
  
  avail_cols <- c(rna_col, prot_col)
  avail_cols <- avail_cols[avail_cols %in% colnames(celltype_scores)]
  
  tmp <- celltype_scores[, avail_cols, drop=FALSE]
  tmp[,] <- lapply(tmp, as.numeric)
  celltype_scores[[comb_col]] <- rowMeans(tmp, na.rm=TRUE)
  
  if (all(is.na(celltype_scores[[comb_col]]))) {
    celltype_scores[[comb_col]] <- NA_real_
  }
  
  ## summary table
  celltype_summary[[ct]] <- data.frame(
    cell_type = ct,
    n_markers_input = length(celltype_markers[[ct]]),
    n_rna_present = res_rna$n,
    n_prot_present = res_prot$n,
    stringsAsFactors = FALSE
  )
  
  ## heatmaps of present markers
  if (res_rna$n >= 5) {
    mat_rna <- X$transcriptome[, res_rna$genes_present, drop=FALSE]
    make_heatmap_nocluster(
      mat_rna,
      group=Y,
      filename_prefix=paste0("HEATMAP_CELLTYPE_", tag, "_RNA_markers"),
      scale_rows=TRUE,
      show_rownames=FALSE,
      cluster_features=TRUE
    )
  }
  
  if (res_prot$n >= 3) {
    mat_prot <- X$proteome[, res_prot$genes_present, drop=FALSE]
    make_heatmap_nocluster(
      mat_prot,
      group=Y,
      filename_prefix=paste0("HEATMAP_CELLTYPE_", tag, "_PROT_markers"),
      scale_rows=TRUE,
      show_rownames=FALSE,
      cluster_features=TRUE
    )
  }
  
  ## boxplot combined score
  df_ct <- celltype_scores %>%
    dplyr::select(sample, genotype, environment, group, all_of(c(rna_col, prot_col, comb_col))) %>%
    dplyr::filter(!is.na(.data[[comb_col]]))
  
  if (nrow(df_ct) > 0) {
    p_ct <- ggplot(df_ct, aes(x=environment, y=.data[[comb_col]], color=environment)) +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(width=0.12, size=1.8) +
      facet_wrap(~ genotype, nrow=1) +
      theme_bw() +
      labs(title=paste0(ct, " | combined marker score"),
           x=NULL,
           y="Combined score (RNA + PROT)") +
      theme(legend.position="none")
    
    save_plot_pdf(paste0("CELLTYPE_", tag, "_combinedScore_byENV_facetGenotype.pdf"), { print(p_ct) },
                  width=12, height=5)
    save_plot_tiff(paste0("CELLTYPE_", tag, "_combinedScore_byENV_facetGenotype.tiff"), { print(p_ct) },
                   width_in=12, height_in=5, dpi=600)
  }
  
  ## ANOVA combined score
  df_aov <- celltype_scores %>%
    dplyr::select(sample, genotype, environment, all_of(comb_col)) %>%
    dplyr::filter(!is.na(.data[[comb_col]]),
                  !is.na(genotype),
                  !is.na(environment))
  
  if (nrow(df_aov) > 0) {
    model_ct <- lm(as.formula(paste0(comb_col, " ~ genotype * environment")), data=df_aov)
    aov_tab  <- anova(model_ct)
    
    terms <- setdiff(rownames(aov_tab), "Residuals")
    if (length(terms) > 0 && "Residuals" %in% rownames(aov_tab)) {
      ss_res <- as.numeric(aov_tab["Residuals","Sum Sq"])
      ss_terms <- as.numeric(aov_tab[terms,"Sum Sq"])
      names(ss_terms) <- terms
      partial_eta2 <- ss_terms / (ss_terms + ss_res)
      
      aov_out <- data.frame(
        cell_type = ct,
        term = terms,
        p_value = as.numeric(aov_tab[terms,"Pr(>F)"]),
        partial_eta2 = as.numeric(partial_eta2),
        stringsAsFactors = FALSE
      )
      
      celltype_anova_list[[ct]] <- aov_out
      
      write_csv(
        aov_out,
        paste0("ANOVA_CELLTYPE_", tag, "_combinedScore_genotype_environment.csv")
      )
    }
  }
}

## save main tables
celltype_summary_df <- do.call(rbind, celltype_summary)
write_csv(celltype_summary_df, "SUMMARY_celltype_marker_detection_RNA_PROT.csv")

write_csv(celltype_scores, "CELLTYPE_scores_all_samples.csv")

if (length(celltype_anova_list) > 0) {
  celltype_anova_df <- do.call(rbind, celltype_anova_list)
  write_csv(celltype_anova_df, "SUMMARY_celltype_ANOVA_combinedScores.csv")
}

## =========================
## 15) OVERVIEW HEATMAPS OF CELL-TYPE SCORES
## =========================
score_cols_comb <- grep("^CT_.*_COMBINED$", colnames(celltype_scores), value=TRUE)
score_cols_rna  <- grep("^CT_.*_RNA$", colnames(celltype_scores), value=TRUE)
score_cols_prot <- grep("^CT_.*_PROT$", colnames(celltype_scores), value=TRUE)

if (length(score_cols_comb) >= 2) {
  mat_comb <- as.matrix(celltype_scores[, score_cols_comb, drop=FALSE])
  rownames(mat_comb) <- celltype_scores$sample
  make_heatmap_nocluster(
    mat_comb,
    group=Y,
    filename_prefix="HEATMAP_allCellTypes_COMBINED_scores",
    scale_rows=TRUE,
    show_rownames=TRUE,
    cluster_features=TRUE
  )
}

if (length(score_cols_rna) >= 2) {
  mat_rna_score <- as.matrix(celltype_scores[, score_cols_rna, drop=FALSE])
  rownames(mat_rna_score) <- celltype_scores$sample
  make_heatmap_nocluster(
    mat_rna_score,
    group=Y,
    filename_prefix="HEATMAP_allCellTypes_RNA_scores",
    scale_rows=TRUE,
    show_rownames=TRUE,
    cluster_features=TRUE
  )
}

if (length(score_cols_prot) >= 2) {
  mat_prot_score <- as.matrix(celltype_scores[, score_cols_prot, drop=FALSE])
  rownames(mat_prot_score) <- celltype_scores$sample
  make_heatmap_nocluster(
    mat_prot_score,
    group=Y,
    filename_prefix="HEATMAP_allCellTypes_PROT_scores",
    scale_rows=TRUE,
    show_rownames=TRUE,
    cluster_features=TRUE
  )
}

cat("\nDONE: cell-type marker analysis completed.\n")
cat("Main outputs:\n")
cat("- CELLTYPE_scores_all_samples.csv\n")
cat("- SUMMARY_celltype_marker_detection_RNA_PROT.csv\n")
cat("- SUMMARY_celltype_ANOVA_combinedScores.csv\n")
cat("- heatmaps and boxplots per cell type\n")

## =========================================================
## 16) DELTA CORRELATION (ISS - GROUND) BY GENOTYPE
##     FINAL CLEAN VERSION FOR MAIN FIGURE
##     - split by genotype
##     - fixed row/column order
##     - DIABLO-like correlation palette (grey20 -> white -> chocolate3)
##     - cleaner labels for publication
## =========================================================

## -------------------------
## A) choose RBPs to display
## -------------------------
rbp_targets <- c("ELAVL1", "RBFOX2", "CIRBP", "FUS", "HNRNPA1", "RBM3", "TARDBP")

## use proteomics as the RBP source
rbp_mat <- X$proteome

find_feature_caseflex <- function(target, cols) {
  if (target %in% cols) return(target)
  cols_norm   <- toupper(gsub("[\\.-]", "_", cols))
  target_norm <- toupper(gsub("[\\.-]", "_", target))
  hit <- cols[cols_norm == target_norm]
  if (length(hit) > 0) return(hit[1])
  hit <- cols[grepl(target_norm, cols_norm, fixed = TRUE)]
  if (length(hit) > 0) return(hit[1])
  return(NA_character_)
}

rbp_found <- sapply(rbp_targets, find_feature_caseflex, cols = colnames(rbp_mat))
rbp_found <- rbp_found[!is.na(rbp_found)]

cat("\nRBPs found in proteomics:\n")
print(rbp_found)

## -------------------------
## B) compute pathway scores
## -------------------------
pathway_scores_df <- data.frame(sample = rownames(X$transcriptome), stringsAsFactors = FALSE)
rownames(pathway_scores_df) <- pathway_scores_df$sample

for (pw in names(pathways_list)) {
  min_g <- dplyr::case_when(
    pw == "Cellular Senescence" ~ 10,
    pw == "METABOLISM_OF_RNA" ~ 15,
    pw == "NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION" ~ 6,
    pw %in% c(
      "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
      "OXIDATIVE_STRESS_RESPONSE",
      "AUTOPHAGY_LYSOSOME",
      "SYNAPTIC_VESICLE_CYCLE"
    ) ~ 4,
    TRUE ~ 6
  )
  
  res_pw <- compute_proxy_score(X$transcriptome, pathways_list[[pw]], min_genes = min_g)
  pathway_scores_df[[pw]] <- res_pw$score
}

## -------------------------
## C) combine metadata + RBPs + pathway scores
## -------------------------
delta_df <- meta
delta_df$sample <- rownames(delta_df)

for (rbp in names(rbp_found)) {
  delta_df[[rbp]] <- rbp_mat[delta_df$sample, rbp_found[[rbp]]]
}

delta_df <- merge(delta_df, pathway_scores_df, by = "sample", all.x = TRUE, sort = FALSE)
delta_df <- delta_df[match(meta$sample, delta_df$sample), ]

## -------------------------
## D) safe Spearman
## -------------------------
safe_spearman <- function(x, y, min_n = 3) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < min_n) return(NA_real_)
  if (length(unique(x)) < 2) return(NA_real_)
  if (length(unique(y)) < 2) return(NA_real_)
  
  suppressWarnings(cor(x, y, method = "spearman"))
}

## -------------------------
## E) compute delta correlation by genotype
## -------------------------
genotypes_to_plot <- c("WT", "APP", "PSEN1", "PSEN2")

pathways_to_plot <- c(
  "NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
  "SYNAPTIC_VESICLE_CYCLE",
  "AUTOPHAGY_LYSOSOME",
  "METABOLISM_OF_RNA",
  "Cellular Senescence",
  "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
  "OXIDATIVE_STRESS_RESPONSE"
)

pretty_pw <- c(
  "NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION" = "Postsynaptic\nsignaling",
  "SYNAPTIC_VESICLE_CYCLE" = "Synaptic vesicle\ncycle",
  "AUTOPHAGY_LYSOSOME" = "Autophagy/\nlysosome",
  "METABOLISM_OF_RNA" = "Metabolism\nof RNA",
  "Cellular Senescence" = "Cellular\nsenescence",
  "TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES" = "TP53 cell\ndeath",
  "OXIDATIVE_STRESS_RESPONSE" = "Oxidative stress\nresponse"
)

delta_long <- list()

for (gt in genotypes_to_plot) {
  df_gt <- delta_df[delta_df$genotype == gt, , drop = FALSE]
  
  for (rbp in names(rbp_found)) {
    for (pw in pathways_to_plot) {
      
      df_ground <- df_gt[df_gt$environment == "GROUND", , drop = FALSE]
      df_iss    <- df_gt[df_gt$environment == "ISS", , drop = FALSE]
      
      rho_ground <- safe_spearman(df_ground[[rbp]], df_ground[[pw]])
      rho_iss    <- safe_spearman(df_iss[[rbp]],    df_iss[[pw]])
      delta_rho  <- rho_iss - rho_ground
      
      delta_long[[length(delta_long) + 1]] <- data.frame(
        genotype   = gt,
        RBP        = rbp,
        pathway    = pw,
        rho_ground = rho_ground,
        rho_iss    = rho_iss,
        delta_rho  = delta_rho,
        stringsAsFactors = FALSE
      )
    }
  }
}

delta_long <- do.call(rbind, delta_long)
write_csv(delta_long, "RBP_pathway_deltaCorrelation_byGenotype.csv")

## -------------------------
## F) plot as faceted heatmap
## -------------------------
plot_df <- delta_long
plot_df$pathway_label <- pretty_pw[plot_df$pathway]

plot_df$RBP <- factor(plot_df$RBP,
                      levels = c("TARDBP", "RBM3", "HNRNPA1", "FUS", "CIRBP", "ELAVL1", "RBFOX2"))

plot_df$pathway_label <- factor(plot_df$pathway_label,
                                levels = c(
                                  "Postsynaptic\nsignaling",
                                  "Synaptic vesicle\ncycle",
                                  "Autophagy/\nlysosome",
                                  "Metabolism\nof RNA",
                                  "Cellular\nsenescence",
                                  "TP53 cell\ndeath",
                                  "Oxidative stress\nresponse"
                                ))

plot_df$genotype <- factor(plot_df$genotype, levels = genotypes_to_plot)

## fixed symmetric scale, cleaner and closer to DIABLO correlation colors
p_delta_by_gt <- ggplot(plot_df, aes(x = pathway_label, y = RBP, fill = delta_rho)) +
  geom_tile(color = NA) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_gradient2(
    low = "#6F97C4",
    mid = "#F3E7A6",
    high = "#E36C4A",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = expression(Delta*rho~"(ISS - Ground)")
  ) +
  labs(
    title = "RBP–pathway delta correlations by genotype",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(face = "italic", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

save_plot_pdf("RBP_pathway_deltaCorrelation_byGenotype_faceted_FINAL.pdf", { print(p_delta_by_gt) },
              width = 14, height = 4.8)

save_plot_tiff("RBP_pathway_deltaCorrelation_byGenotype_faceted_FINAL.tiff", { print(p_delta_by_gt) },
               width_in = 14, height_in = 4.8, dpi = 600)

# PATCH NOTES FOR YOUR FULL SCRIPT

# 1) In the circos section, replace:
# cutoff_circos <- 0.7
# with:
cutoff_circos <- 0.75

# 2) In BOTH circosPlot calls, replace:
# size.labels=1.5
# with:
size.labels=1.1

# 3) Optional cleaner theme for CIRBP module plot (section 11, p6B):
# replace the p6B block with:

p6B <- ggplot(scores, aes(x=environment, y=CIRBP_multiomic_module, color=environment)) +
  geom_boxplot(outlier.shape=NA, width=0.55) +
  geom_jitter(width=0.10, size=1.8) +
  facet_wrap(~ genotype, nrow=1) +
  scale_color_manual(values=c("GROUND"="#F8766D","ISS"="#00BFC4")) +
  labs(title="CIRBP multi-omic module by environment", x=NULL, y="Module score (z)") +
  theme_bw(base_size=11) +
  theme(
    strip.background = element_rect(fill="grey95", color="grey70"),
    strip.text = element_text(face="bold"),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(face="bold", hjust=0.5)
  )

# 4) Optional cleaner theme for any RBM3 module plot:
# use the same theme structure as above, just change title/y variable.

# 5) Replace your current ENTIRE section 16 with the code in the companion file:
# "section16_delta_correlation_final.R"

## =========================
## 13) REPORT
## =========================
report_path <- file.path(outdir, "DIABLO_report_README.txt")
report_lines <- c(
  "mixOmics DIABLO — PLAN B — CLEAN UPDATED (manual pathways only)",
  "",
  "Preprocessing:",
  "- drop zero-variance features",
  "- median imputation per feature",
  "- scale features",
  "- near-zero variance removal",
  "",
  "Model:",
  "- block.splsda with ncomp=2, keepX fixed, design=0.5 between blocks",
  "",
  "Figures:",
  "- DIABLO individuals/variables/block cor/arrow/CIM/circos",
  "- heatmaps (no sample clustering)",
  "- env-focused heatmaps (genotype removed; env preserved)",
  "- Fig6B CIRBP module by ENV (facet genotype)",
  "- Fig6C CIRBP vs manual pathway proxy scores (scatter)",
  "- ExtFig6-3 CIRBP-linked module heatmaps",
  "- Manual pathway analyses (ANOVA/posthoc/plots/corr) for the 5 sets",
  "",
  "Outputs:",
  paste0("- PDF figures:  ", dir_fig_pdf),
  paste0("- TIFF figures: ", dir_fig_tiff),
  paste0("- CSV tables:   ", dir_tables),
  paste0("- Excel workbook: ", xlsx_path),
  paste0("- Saved model: ", file.path(outdir, "DIABLO_model_planB.rds")),
  ""
)
writeLines(report_lines, report_path)

cat("\nDONE ✅\n")
cat("PDF figures:  ", dir_fig_pdf, "\n", sep="")
cat("TIFF figures: ", dir_fig_tiff, "\n", sep="")
cat("CSV tables:   ", dir_tables, "\n", sep="")
cat("Excel workbook:\n", xlsx_path, "\n", sep="")
cat("Report:\n", report_path, "\n", sep="")