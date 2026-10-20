# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

## `DESeq2` differential expression analysis configuration

To successfully run the differential expression analysis, you will need to tell DESeq2 which sample annotations to use (annotations are columns in the `samples.tsv` file described below).
This is done in the `config.yaml` file with the entries under `diffexp:`.
The comments for the entries should give all the necessary infos and linkouts.
But if in doubt, please also consult the [`DESeq2` manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

## Starting from a precomputed count matrix

Set `counts:` in `config/config.yaml` to a tab-separated gene count matrix (first column `gene`, one column per `sample_name` of the sample sheet) to skip trimming, alignment and read-level QC.
Columns that are not in the sample sheet are ignored, so one matrix can serve several sample sheets.
Leave `counts:` empty to run the workflow from the FASTQ files in the unit sheet.

## Interaction contrasts

Besides the `variable_of_interest`/`level_of_interest` form, a contrast can be given as a string holding a DESeq2 `contrast` expression, e.g. `'list(c("ENVIRONMENTISS.DISEASEAD"))'` to test an interaction coefficient, or `'list(c("ENVIRONMENT_ISS_vs_GROUND", "ENVIRONMENTISS.DISEASEAD"))'` for the effect of one variable within a non-reference level of another.
Coefficient names are those returned by `DESeq2::resultsNames()` for the fitted model.

# Sample and unit setup

The sample and unit setup is specified via tab-separated tabular files (`.tsv`).
Missing values can be specified by empty columns or by writing `NA`.

## sample sheet

The default sample sheet is `config/samples.tsv` (as configured in `config/config.yaml`).
Each sample refers to an actual physical sample, and replicates (both biological and technical) may be specified as separate samples.
For each sample, you will always have to specify a `sample_name`.
In addition, all `variables_of_interest` and `batch_effects` specified in the `config/config.yaml` under the `diffexp:` entry, will have to have corresponding columns in the `config/samples.tsv`.
Finally, the sample sheet can contain any number of additional columns.
So if in doubt about whether you might at some point need some metadata you already have at hand, just put it into the sample sheet already---your future self will thank you.

## unit sheet

The default unit sheet is `config/units.tsv` (as configured in `config/config.yaml`).
For each sample, add one or more sequencing units (for example if you have several runs or lanes per sample).

### `.fastq` file source

For each unit, you will have to define a source for your `.fastq` files.
This can be done via the columns `fq1`, `fq2` and `sra`, with either of:
1. A single `.fastq` file for single-end reads (`fq1` column only; `fq2` and `sra` columns present, but empty).
  The entry can be any path on your system, but we suggest something like a `raw/` data directory within your analysis directory.
2. Two `.fastq` files for paired-end reads (columns `fq1` and `fq2`; column `sra` present, but empty).
  As for the `fq1` column, the `fq2` column can also point to anywhere on your system.
3. A sequence read archive (SRA) accession number (`sra` column only; `fq1` and `fq2` columns present, but empty).
  The workflow will automatically download the corresponding `.fastq` data (currently assumed to be paired-end).
  The accession numbers usually start with SRR or ERR and you can find accession numbers for studies of interest with the [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/).
If both local files and an SRA accession are specified for the same unit, the local files will be used.

### adapter trimming

If you set `trimming: activate:` in the `config/config.yaml` to `True`, reads are trimmed with [Trim Galore](https://github.com/FelixKrueger/TrimGalore) using the Illumina adapter and a quality cutoff of 20 (`--illumina -q 20`) before alignment.

### strandedness of library preparation protocol

To get the correct `geneCounts` from `STAR` output, you can provide information on the strandedness of the library preparation protocol used for a unit.
`STAR` can produce counts for unstranded (`none` - this is the default), forward oriented (`yes`) and reverse oriented (`reverse`) protocols.  
Enter the respective value into a `strandedness` column in the `units.tsv` file.
