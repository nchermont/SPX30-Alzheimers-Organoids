# =========================
# REACTOME ENRICHMENT PER CLUSTER
# from heatmap_data.csv
# =========================

# =========================
# 1. PACKAGES
# =========================

if(!require("BiocManager")) install.packages("BiocManager")

bioc_packages <- c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "enrichplot"
)

for(pkg in bioc_packages){
  if(!require(pkg, character.only = TRUE)){
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    library(pkg, character.only = TRUE)
  }
}

cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "openxlsx",
  "stringr"
)

for(pkg in cran_packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =========================
# 2. FILES
# =========================

heatmap_file <- ""

output_dir <- ""
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Saving files in:\n", normalizePath(output_dir), "\n")

# =========================
# 3. LOAD DATA
# =========================

heatmap <- read.csv(
  heatmap_file,
  check.names = FALSE
)

colnames(heatmap) <- trimws(colnames(heatmap))
heatmap <- heatmap[, colnames(heatmap) != ""]

# =========================
# 4. PREPARE GENE LISTS
# using Ensembl IDs
# =========================

cluster_gene_df <- heatmap %>%
  dplyr::filter(!is.na(ensembl_ID), ensembl_ID != "") %>%
  dplyr::select(ensembl_ID, symbol, cluster) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    ensembl_ID_clean = gsub("\\..*", "", ensembl_ID),
    Cluster = paste0("Cluster ", cluster)
  )

# =========================
# 5. ENSEMBL TO ENTREZ ID
# =========================

gene_mapping <- bitr(
  cluster_gene_df$ensembl_ID_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cluster_gene_mapped <- cluster_gene_df %>%
  dplyr::inner_join(
    gene_mapping,
    by = c("ensembl_ID_clean" = "ENSEMBL")
  )

cat("Mapped genes:\n")
print(nrow(cluster_gene_mapped))

# =========================
# 6. REACTOME ENRICHMENT PER CLUSTER
# =========================

reactome_results <- list()
clusters <- unique(cluster_gene_mapped$Cluster)

for(cl in clusters){
  
  message("Running Reactome enrichment for ", cl)
  
  entrez_genes <- cluster_gene_mapped %>%
    dplyr::filter(Cluster == cl) %>%
    dplyr::pull(ENTREZID) %>%
    unique()
  
  enrich_res <- enrichPathway(
    gene = entrez_genes,
    organism = "human",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  reactome_results[[cl]] <- enrich_res
  result_df <- as.data.frame(enrich_res)
  
  write.xlsx(
    result_df,
    file.path(output_dir, paste0(gsub(" ", "_", cl), "_Reactome_results.xlsx")),
    rowNames = FALSE
  )
  
  if(nrow(result_df) > 0){
    
    result_top6 <- result_df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = 6) %>%
      dplyr::mutate(
        Description_wrapped = stringr::str_wrap(Description, width = 42),
        logFDR = -log10(p.adjust),
        GeneRatio_decimal = as.numeric(sub("/.*", "", GeneRatio)) /
          as.numeric(sub(".*/", "", GeneRatio))
      )
    
    result_top6$Description_wrapped <- factor(
      result_top6$Description_wrapped,
      levels = rev(result_top6$Description_wrapped)
    )
    
    # =========================
    # DOTPLOT PER CLUSTER
    # =========================
    
    p_dot <- ggplot(
      result_top6,
      aes(
        x = GeneRatio_decimal,
        y = Description_wrapped,
        size = Count,
        color = logFDR
      )
    ) +
      geom_point(alpha = 0.9) +
      scale_color_gradient(
        low = "#B7E4C7",
        high = "#1B4332"
      ) +
      scale_size(range = c(2.5, 7)) +
      scale_y_discrete(expand = c(0.02, 0.02)) +
      theme_classic(base_size = 13) +
      theme(
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold", lineheight = 0.9),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 17)
      ) +
      labs(
        title = paste0(cl, " - Reactome enrichment"),
        x = "GeneRatio",
        y = "",
        size = "Gene count",
        color = "-log10(FDR)"
      )
    
    ggsave(
      file.path(output_dir, paste0(gsub(" ", "_", cl), "_Reactome_dotplot.pdf")),
      p_dot,
      width = 10,
      height = 4
    )
    
    ggsave(
      file.path(output_dir, paste0(gsub(" ", "_", cl), "_Reactome_dotplot.png")),
      p_dot,
      width = 10,
      height = 4,
      dpi = 300
    )
    
    # =========================
    # BARPLOT PER CLUSTER
    # =========================
    
    p_bar <- ggplot(
      result_top6,
      aes(
        x = logFDR,
        y = Description_wrapped,
        fill = logFDR
      )
    ) +
      geom_col(width = 0.65) +
      scale_fill_gradient(
        low = "#B7E4C7",
        high = "#1B4332"
      ) +
      scale_y_discrete(expand = c(0.02, 0.02)) +
      theme_classic(base_size = 13) +
      theme(
        axis.text.x = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold", lineheight = 0.9),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 17)
      ) +
      labs(
        title = paste0(cl, " - Reactome enrichment"),
        x = "-log10(FDR)",
        y = "",
        fill = "-log10(FDR)"
      )
    
    ggsave(
      file.path(output_dir, paste0(gsub(" ", "_", cl), "_Reactome_barplot.pdf")),
      p_bar,
      width = 10,
      height = 4
    )
    
    ggsave(
      file.path(output_dir, paste0(gsub(" ", "_", cl), "_Reactome_barplot.png")),
      p_bar,
      width = 10,
      height = 4,
      dpi = 300
    )
  }
}

# =========================
# 7. COMBINE ALL RESULTS
# =========================

all_reactome_df <- bind_rows(
  lapply(names(reactome_results), function(cl){
    df <- as.data.frame(reactome_results[[cl]])
    if(nrow(df) > 0){
      df$Cluster <- cl
      return(df)
    }
  })
)

write.xlsx(
  all_reactome_df,
  file.path(output_dir, "All_clusters_Reactome_enrichment_combined.xlsx"),
  rowNames = FALSE
)

# =========================
# 8. COMBINED DOTPLOT
# =========================

if(nrow(all_reactome_df) > 0){
  
  top_terms <- all_reactome_df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 6) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Description_wrapped = stringr::str_wrap(Description, width = 36),
      logFDR = -log10(p.adjust),
      GeneRatio_decimal = as.numeric(sub("/.*", "", GeneRatio)) /
        as.numeric(sub(".*/", "", GeneRatio))
    )
  
  p_combined <- ggplot(
    top_terms,
    aes(
      x = Cluster,
      y = reorder(Description_wrapped, -p.adjust),
      size = Count,
      color = logFDR
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#B7E4C7",
      high = "#1B4332"
    ) +
    scale_size(range = c(2.5, 7)) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 13,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 10,
        face = "bold",
        lineheight = 0.9
      ),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
    ) +
    labs(
      title = "Top Reactome pathways across expression clusters",
      x = "",
      y = "",
      size = "Gene count",
      color = "-log10(FDR)"
    )
  
  print(p_combined)
  
  ggsave(
    file.path(output_dir, "All_clusters_Reactome_combined_dotplot.pdf"),
    p_combined,
    width = 7,
    height = 15
  )
  
  ggsave(
    file.path(output_dir, "All_clusters_Reactome_combined_dotplot.png"),
    p_combined,
    width = 7,
    height = 15,
    dpi = 300
  )
}

# =========================
# 9. DONE
# =========================

cat("\nNumber of enriched terms found:\n")
print(table(all_reactome_df$Cluster))

system(paste("open", shQuote(output_dir)))
