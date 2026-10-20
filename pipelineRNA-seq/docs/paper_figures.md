# Paper tables and figures

Rule `paper_figures` (`workflow/rules/paper.smk`, script
`workflow/scripts/paper_figures.R`) builds every table and figure of the paper
from the gene count matrix. It is enabled by the `paper:` section of
`config/config.yaml` and runs with the rest of the workflow:

```
snakemake --cores all --use-conda
```

Outputs are written to `results/paper_figures/`:

| Path | Content |
|------|---------|
| `pca_plot*.png` | PCA of all samples and of each AD line with WTAD |
| `heatmap_top1000_variable_genes.*` | 1000 most variable genes |
| `interaction_<G>_vs_WTAD_ground_vs_iss.csv`, `volcano_interaction_*`, `heatmap_interaction_*`, `heatmap_app.pdf` | genotype × environment interaction, `<G>` = APP, PSEN1, PSEN2 |
| `<G>_vs_WTAD_ISS_enrichGO_BP.csv`, `*_dotplot_BP.png`, `*_barplot_BP.png` | GO biological-process enrichment of the interaction genes |
| `APP/`, `PSEN1/`, `PSEN2/`, `control/` | ISS vs GROUND within each line (`control` = WTAD) |
| `AD-WTAD/` | ISS vs GROUND in the pooled AD lines and in WTAD, gene lists and Venn diagrams |
| `ALL/` | ISS vs GROUND and AD vs WTAD over all samples, gene lists and Venn diagrams |
| `DEG_summary_individual_compare.csv` | number of tested genes and DEGs per ISS vs GROUND comparison |
| `normalized_counts_*.csv`, `vst_transformed_counts_*.csv` | DESeq2-normalised counts and row z-scores of VST values |
| `sample_names.tsv`, `sessionInfo.txt` | sample-name schemes and software versions |

## Inputs

| File | Content | Provenance |
|------|---------|------------|
| `data/counts/all.tsv` | STAR gene counts, 46 samples (not in the repository) | this workflow, run separately on the AD and wtAD samples and merged with `scripts/merge_inputs.py` |
| `config/samples.tsv`, `config/units.tsv` | sample and sequencing-unit sheets | `scripts/merge_inputs.py` |
| `data/gene_annotation.tsv` | gene names, frozen | `scripts/freeze_gene_annotation.py` |

The sequencing data are not part of the repository. `results/counts/all.tsv`
is copied from `data/counts/all.tsv` (`counts:` in `config/config.yaml`),
which has to be provided: a tab-separated matrix with a `gene` column of
Ensembl ids followed by one column per `sample_name` of `config/samples.tsv`.
With `counts: ""` the counts are instead computed from the FASTQ files listed
in `config/units.tsv`, expected in `data/fastq/` (paired-end,
`<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`).

Gene names were retrieved in August 2025 from biomaRt (Ensembl
`external_gene_name`; interaction, normalised-count and `ALL/` tables) and
from org.Hs.eg.db (per-line ISS vs GROUND tables). They are stored in
`data/gene_annotation.tsv` so that the tables can be rebuilt offline with the
same names. Genes without a name are reported by their Ensembl id in the
tables and are not labelled in the volcano plots.

## Methods

Samples: APP, PSEN1, PSEN2 and WTAD lines grown on GROUND and on the ISS;
3 biological replicates (2 for GROUND PSEN1), each sequenced on two lanes
(L005, L006) that enter the analysis as separate samples (46 in total).

All tests are DESeq2 Wald tests without log2 fold-change shrinkage, with
default independent filtering (alpha = 0.1). Differentially expressed genes
(DEGs): padj < 0.05 and |log2 fold change| > 1.

| Output | Samples | Genes kept | Model / transformation | Result |
|--------|---------|------------|------------------------|--------|
| `interaction_<G>_...csv`, interaction heatmaps and volcano plots, GO | all 46 | row sum > 10 | `~condition * environment` (condition: WTAD, APP, PSEN1, PSEN2) | coefficient `condition<G>.environmentISS`: ISS effect in `<G>` minus ISS effect in WTAD |
| `APP/`, `PSEN1/`, `PSEN2/`, `control/`, `AD-WTAD/*_iss_vs_ground_all.csv`, `DEG_summary` | GROUND and ISS samples of the line(s) | genes with any count | `~environment` | ISS vs GROUND; padj of genes removed by independent filtering written as 1 |
| `ALL/DESeq2_ISS_vs_GROUND_results_ordered.csv` | all 46 | row sum > 10 | `~condition * environment` | `environment_ISS_vs_GROUND`: the ISS effect in the reference line, WTAD |
| `ALL/DESeq2_AD_vs_WTAD_results_ordered.csv` | all 46 | row sum > 10 | `~disease` (AD, WTAD) | AD vs WTAD |
| `pca_plot.png`, `heatmap_top1000_*` | all 46 | row sum > 10 | `vst(blind = TRUE)`; top 500 (PCA) or 1000 (heatmap) genes by variance; heatmap rows z-scored | |
| `pca_plot_<G>.png` | `<G>` and WTAD | row sum > 10 | `vst(blind = TRUE)` of these samples, top 500 genes | |
| `heatmap_interaction_*`, `heatmap_app.pdf` | `<G>` and WTAD shown | row sum > 10 | `vst(blind = FALSE)` of the interaction model, centred on the row mean of all 46 samples | 30 interaction genes with the lowest padj |
| `*_enrichGO_BP.csv` | | universe: all genes of the interaction table | clusterProfiler `enrichGO(ont = "BP", keyType = "SYMBOL")`, BH adjustment, no cut-off | interaction DEGs of `<G>` |

Volcano plots of single lines are capped at |log2FC| = 5 and −log10 FDR = 50
and label the 10 most significant up- and down-regulated genes that have a
symbol. Interaction volcano plots label the genes marked in the highlighted
heatmaps.

## Sample names

`results/paper_figures/sample_names.tsv` maps three naming schemes:

| Scheme | Example | Used in |
|--------|---------|---------|
| sample sheet | `GROUND_PSEN1_B_S37_L005` | all tables except below |
| numbered | `GROUND_PSEN1_3` | `vst_transformed_counts_all.csv`, `heatmap_app.pdf` |
| standard | `GROUND_PSEN1_B_L005` | environment, line, replicate, lane |

Numbered names count replicates A, B, C × lanes L005, L006 in that order
within each environment and line.

## Colours

| Category | Colour |
|----------|--------|
| upregulated (higher in ISS, or in AD) | blue `#0000FF` (light `#7F7FFF` for a second "up" set in Venn diagrams) |
| downregulated | orange `#FFA500` (light `#FFD27F` for a second "down" set) |
| not significant | grey |
| heatmap scale | orange (low) – white – blue (high) |
| GO plots, p.adjust | orange (lowest) – blue |
| WTAD / APP / PSEN1 / PSEN2 | grey `#7F7F7F` / orange `#E69F00` / green `#009E73` / pink `#CC79A7` |
| GROUND / ISS | brown `#8C6D31` / sky blue `#56B4E9` |
