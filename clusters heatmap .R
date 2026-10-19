# =========================
# CLUSTER LINE PLOT
# Mean of triplicates
# One plot per cluster
# =========================

packages <- c("tidyverse", "ggplot2", "patchwork")

for(pkg in packages){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =========================
# FILES
# =========================

heatmap_file <- "/Users/namoreira/Library/CloudStorage/GoogleDrive-nchmoreira0@gmail.com/My Drive/POSDOC UCSD/Experimentos/space project/SPX31 data/BULK PROTEOMIC SPX31/HEATMAP TOP 1000 GERAL/heatmap_data (5).csv"

output_dir <- "Cluster_line_plots_mean_triplicatesall
"
dir.create(output_dir, showWarnings = FALSE)

# =========================
# LOAD DATA
# =========================

heatmap <- read.csv(
  heatmap_file,
  check.names = FALSE
)

colnames(heatmap) <- trimws(colnames(heatmap))

# Remove empty/unnamed columns
heatmap <- heatmap[, colnames(heatmap) != ""]

print(colnames(heatmap))

# =========================
# COLUMNS
# =========================

gene_col <- "symbol"
cluster_col <- "cluster"

sample_cols <- grep(
  "Ground|ISS",
  colnames(heatmap),
  value = TRUE
)

print(sample_cols)

# =========================
# LONG FORMAT
# =========================

long_df <- heatmap %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  mutate(
    Expression = as.numeric(Expression),
    Cluster = paste0("Cluster ", .data[[cluster_col]]),
    Gene = .data[[gene_col]]
  )

# =========================
# CREATE BIOLOGICAL GROUPS
# =========================

long_df <- long_df %>%
  mutate(
    Group = case_when(
      grepl("NDC2_Ground", Sample) ~ "NDC2_Ground",
      grepl("NDC2_ISS", Sample) ~ "NDC2_ISS",
      grepl("sAD1_Ground", Sample) ~ "sAD1_Ground",
      grepl("sAD1_ISS", Sample) ~ "sAD1_ISS",
      grepl("sAD2_Ground", Sample) ~ "sAD2_Ground",
      grepl("sAD2_ISS", Sample) ~ "sAD2_ISS",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# =========================
# MEAN OF TRIPLICATES
# =========================

mean_triplicate_df <- long_df %>%
  group_by(Gene, Cluster, Group) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# GROUP ORDER
# =========================

group_order <- c(
  "NDC2_Ground",
  "NDC2_ISS",
  "sAD1_Ground",
  "sAD1_ISS",
  "sAD2_Ground",
  "sAD2_ISS"
)

mean_triplicate_df <- mean_triplicate_df %>%
  mutate(
    Group = factor(Group, levels = group_order),
    Group_index = as.numeric(Group)
  )

# =========================
# Z-SCORE BY GENE
# =========================

mean_triplicate_df <- mean_triplicate_df %>%
  group_by(Gene) %>%
  mutate(
    Zscore = as.numeric(scale(Mean_Expression))
  ) %>%
  ungroup()

# =========================
# MEAN LINE PER CLUSTER
# =========================

cluster_mean_df <- mean_triplicate_df %>%
  group_by(Cluster, Group, Group_index) %>%
  summarise(
    Mean_Z = mean(Zscore, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# PLOT ALL CLUSTERS TOGETHER
# =========================

p_all <- ggplot() +
  geom_line(
    data = mean_triplicate_df,
    aes(
      x = Group_index,
      y = Zscore,
      group = Gene
    ),
    color = "#F4A300",
    alpha = 0.22,
    linewidth = 0.35
  ) +
  geom_line(
    data = cluster_mean_df,
    aes(
      x = Group_index,
      y = Mean_Z,
      group = Cluster
    ),
    color = "#B22222",
    linewidth = 1.6
  ) +
  geom_point(
    data = cluster_mean_df,
    aes(
      x = Group_index,
      y = Mean_Z
    ),
    color = "#B22222",
    size = 2.3
  ) +
  facet_wrap(~Cluster, ncol = 2) +
  scale_x_continuous(
    breaks = 1:length(group_order),
    labels = group_order
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 16
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    )
  ) +
  labs(
    title = "Expression dynamics across heatmap clusters",
    x = "",
    y = "Mean triplicate Z-score"
  )

print(p_all)

ggsave(
  file.path(output_dir, "All_clusters_mean_triplicates.pdf"),
  p_all,
  width = 12,
  height = 8
)

ggsave(
  file.path(output_dir, "All_clusters_mean_triplicates.png"),
  p_all,
  width = 12,
  height = 8,
  dpi = 300
)

# =========================
# SAVE ONE PLOT PER CLUSTER
# =========================

clusters <- unique(mean_triplicate_df$Cluster)

for(cl in clusters){
  
  plot_data <- mean_triplicate_df %>%
    filter(Cluster == cl)
  
  mean_data <- cluster_mean_df %>%
    filter(Cluster == cl)
  
  p_cluster <- ggplot() +
    geom_line(
      data = plot_data,
      aes(
        x = Group_index,
        y = Zscore,
        group = Gene
      ),
      color = "#F4A300",
      alpha = 0.25,
      linewidth = 0.35
    ) +
    geom_line(
      data = mean_data,
      aes(
        x = Group_index,
        y = Mean_Z,
        group = 1
      ),
      color = "#B22222",
      linewidth = 1.8
    ) +
    geom_point(
      data = mean_data,
      aes(
        x = Group_index,
        y = Mean_Z
      ),
      color = "#B22222",
      size = 2.8
    ) +
    scale_x_continuous(
      breaks = 1:length(group_order),
      labels = group_order
    ) +
    theme_classic(base_size = 14) +
    theme(
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 15,
          face = "bold"
        ),
        
        axis.text.y = element_text(
          size = 15,
          face = "bold"
        ),
        
        axis.title = element_text(
          size = 16,
          face = "bold"
        ),
        
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 21
        )
      )
    ) +
    labs(
      title = cl,
      x = "",
      y = "Mean triplicate Z-score"
    )
  
  print(p_cluster)
  
  file_name <- gsub(" ", "_", cl)
  
  ggsave(
    file.path(output_dir, paste0(file_name, "_mean_triplicates.pdf")),
    p_cluster,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(output_dir, paste0(file_name, "_mean_triplicates.png")),
    p_cluster,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# =========================
# SAVE DATA
# =========================

write.csv(
  long_df,
  file.path(output_dir, "Original_long_data_with_replicates.csv"),
  row.names = FALSE
)

write.csv(
  mean_triplicate_df,
  file.path(output_dir, "Mean_triplicate_gene_cluster_data.csv"),
  row.names = FALSE
)

write.csv(
  cluster_mean_df,
  file.path(output_dir, "Mean_triplicate_cluster_average_line.csv"),
  row.names = FALSE
)

cat("\nAnalysis completed successfully!\n")
cat("Files saved in:", output_dir, "\n")