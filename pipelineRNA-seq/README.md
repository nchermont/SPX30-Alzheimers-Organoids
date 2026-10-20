# Snakemake workflow: pipelineRNA-seq

Bulk RNA-seq workflow for read trimming, alignment, gene counting, quality control
and differential expression analysis, part of
[SPX30-Alzheimers-Organoids](https://github.com/nchermont/SPX30-Alzheimers-Organoids).
It is derived from
[snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
(MIT licence, see `LICENSE`).

| Step | Tool |
|------|------|
| Adapter / quality trimming | [Trim Galore](https://github.com/FelixKrueger/TrimGalore) (`--illumina -q 20`) |
| Alignment and gene counts | [STAR](https://github.com/alexdobin/STAR) (`--quantMode GeneCounts`) against the Ensembl genome and annotation set in `config.yaml` |
| Read-level QC | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [RSeQC](https://rseqc.sourceforge.net/), [MultiQC](https://multiqc.info/) |
| Normalisation | TPM and CPM tables |
| Differential expression | [DESeq2](https://bioconductor.org/packages/DESeq2) (Wald test, `ashr` log2 fold-change shrinkage) |
| Visualisation | PCA, heatmaps, MA and volcano plots |
| Enrichment (optional) | GSEA / Fisher tests on MSigDB and Reactome gene sets |

All software is deployed through the conda environments in `workflow/envs/` and the
pinned [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io) (`v7.0.0`).

## Reproducing the paper analysis

This folder is configured for the paper's RNA-seq data: Alzheimer's
disease lines (APP, PSEN1, PSEN2) and a wild-type line (WTAD) grown on the
ground (GROUND) and on the International Space Station (ISS), 46 samples.

```
git clone https://github.com/nchermont/SPX30-Alzheimers-Organoids.git
cd SPX30-Alzheimers-Organoids/pipelineRNA-seq
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --cores all --use-conda
```

All paths are relative to the `pipelineRNA-seq` folder.

| Path | Content |
|------|---------|
| `config/config.yaml`, `config/samples.tsv`, `config/units.tsv` | paper configuration and sample sheets |
| `data/gene_annotation.tsv` | frozen gene names used in the paper tables |
| `results/paper_figures/` | tables and figures of the paper |
| `results/diffexp/` | environment × disease interaction (`~ENVIRONMENT + DISEASE + ENVIRONMENT:DISEASE`) |
| `scripts/` | how `config/` and `data/` were produced from the per-project runs |

The sequencing data are not part of this repository. Before running, provide
one of:

- the gene count matrix of the 46 samples as `data/counts/all.tsv` (the
  default, `counts: data/counts/all.tsv` in `config/config.yaml`), or
- the FASTQ files listed in `config/units.tsv` in `data/fastq/`, with
  `counts: ""` in `config/config.yaml` to compute the counts from the reads.

Methods, sample names and colours of the paper figures are described in
[`docs/paper_figures.md`](docs/paper_figures.md).

The rest of this README describes how to use the workflow on other data.

## Step 1: Install Snakemake

```
mamba create -c conda-forge -c bioconda --name snakemake snakemake
conda activate snakemake
```

## Step 2: Use the workflow as a module

Create a project directory with a `config/` folder (copy `config/` from this
folder as a template) and a `workflow/Snakefile` that imports this workflow as
a Snakemake module:

```python
configfile: "config/config.yaml"


module pipelineRNA_seq:
    snakefile:
        github("nchermont/SPX30-Alzheimers-Organoids", path="pipelineRNA-seq/workflow/Snakefile", branch="main")
    config:
        config


use rule * from pipelineRNA_seq
```

For a reproducible analysis, replace `branch="main"` by `tag="<release>"` or
`commit="<sha>"`.

## Step 3: Configure the workflow

Edit `config/config.yaml` following the comments in the file. Details on every
option are in [`config/README.md`](config/README.md). The main entries are:

- `samples` – sample sheet (`sample_name` plus one column per variable of interest or batch effect).
- `units` – sequencing units (lanes/runs) per sample, with paths to FASTQ files or SRA accessions.
- `counts` – *optional* path to a precomputed gene count matrix (see below).
- `diffexp` – DESeq2 design (`variables_of_interest`, `batch_effects`, `model`) and `contrasts`.

### Starting from a count matrix

If `counts` is set to a TSV file (first column `gene`, one column per `sample_name`),
trimming, alignment and read-level QC are skipped and the analysis starts from that
matrix. This makes it possible to reproduce the statistical analysis from the
deposited count tables without the raw reads. Leave `counts` empty (or omit it) to
run the full workflow from FASTQ files.

### Contrasts and interaction terms

A contrast is either a simple comparison between two levels of one variable:

```yaml
contrasts:
  ISS-vs-GROUND:
    variable_of_interest: ENVIRONMENT
    level_of_interest: ISS
```

or a DESeq2 `contrast` expression written as a string, which allows testing
interaction terms of a model such as `~ENVIRONMENT * DISEASE`:

```yaml
contrasts:
  interaction_ISS_DISEASE: 'list(c("ENVIRONMENTISS.DISEASEAD"))'
```

The coefficient names follow `DESeq2::resultsNames()`; they are written to the
log `logs/deseq2/init.log`.

## Step 4: Run the workflow

```
snakemake --cores all --use-conda
```

Results are written to `results/`:

| Path | Content |
|------|---------|
| `results/counts/all.tsv` | raw gene counts (technical replicates of a sample summed) |
| `results/counts/normalized_{tpm,cpm}.tsv` | TPM / CPM tables |
| `results/deseq2/normcounts.tsv` | DESeq2 size-factor normalised counts |
| `results/diffexp/{contrast}.diffexp.tsv` | DESeq2 results per contrast |
| `results/diffexp/{contrast}.expressiontable.tsv` | results plus group means and normalised counts |
| `results/pca.*.svg`, `results/heatmap.*.svg` | exploratory plots |
| `results/qc/multiqc_report.html` | read-level QC report |

Files with the `.symbol.tsv` suffix carry an additional gene symbol column.

A self-contained HTML report can be generated with `snakemake --report report.zip`.
