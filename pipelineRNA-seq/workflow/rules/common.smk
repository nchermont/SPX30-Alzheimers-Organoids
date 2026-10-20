import glob

import pandas as pd
#from snakemake.remote import FTP
from snakemake.utils import validate
import subprocess as sp


#ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")

def use_precomputed_counts():
    # a precomputed count matrix replaces trimming, alignment and read-level QC
    return bool(config.get("counts"))


def get_final_output():
    final_output=list()
    if not use_precomputed_counts():
        final_output.append(expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/Aligned.sortedByCoord.out.bam.bai",
            unit=units.itertuples(),
        ))
    if config["options"]["diffexp"]:
        final_output.append(expand(
            "results/diffexp/{contrast}.diffexp.symbol.tsv",
            contrast=config["diffexp"]["contrasts"],
        ))
        final_output.append(expand(
            "results/diffexp/{contrast}.expressiontable.tsv",
            contrast=config["diffexp"]["contrasts"],
        ))
        final_output.append(expand(
            "results/diffexp/{contrast}.expressiontable.symbol.tsv",
            contrast=config["diffexp"]["contrasts"],
        ))
        final_output.append(expand(
            "results/diffexp/{contrast}.volcanoplot.svg",
            contrast=config["diffexp"]["contrasts"],
        ))
        if config["pca"]["activate"]:
            # get all the variables to plot a PCA for
            pca_variables = list(config["diffexp"]["variables_of_interest"])
            if config["diffexp"]["batch_effects"]:
                pca_variables.extend(config["diffexp"]["batch_effects"])
            if config["pca"]["labels"]:
                pca_variables.extend(config["pca"]["labels"])
            final_output.extend(
                expand("results/pca.{variable}.svg", variable=pca_variables)
            )
        if config["heatmap"]["activate"]:
            # get all the variables to plot a PCA for
            heatmap_variables = list(config["diffexp"]["variables_of_interest"])
            if config["diffexp"]["batch_effects"]:
                heatmap_variables.extend(config["diffexp"]["batch_effects"])
            if config["heatmap"]["labels"]:
                heatmap_variables.extend(config["heatmap"]["labels"])
            final_output.extend(
                expand("results/heatmap.{variable}_top20.svg", variable=heatmap_variables)
            )
        final_output.append("results/deseq2/normcounts.symbol.tsv")
        final_output.append("results/counts/all.symbol.tsv")
    if config["options"]["gsea"]:
        final_output.append(expand(
            "results/diffexp/{contrast}.gseares.RDS",
            contrast=config["diffexp"]["contrasts"],
        ))
    if config["options"]["normalization"]:
        final_output.append("results/counts/normalized_tpm.tsv")
    if config.get("paper", {}).get("activate", False):
        final_output.append("results/paper_figures")
    if config["options"]["qc"] and not use_precomputed_counts():
        final_output.append("results/qc/multiqc_report.html")
    return final_output





wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq.gz", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def get_fq(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/trimmed/{sample}_{unit}_{group}.fastq.gz",
                        group=["R1", "R2"],
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {"fq1": "results/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)}
    else:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit)]
        if pd.isna(u["fq1"]):
            # SRA sample (always paired-end for now)
            accession = u["sra"]
            if sra_is_paired_end(accession):
                return dict(
                    zip(
                        ["fq1", "fq2"],
                        expand(
                            "sra/pe/{accession}_{group}.fastq.gz",
                            accession=accession,
                            group=["1", "2"],
                        ),
                    )
                )
            else:
                return {"fq1": f"sra/se/{accession}.fastq.gz"}
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else:
            return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}

def get_fq_trimgalore(wildcards):
    # no trimming, use raw reads
    u = units.loc[(wildcards.sample, wildcards.unit)]
    if pd.isna(u["fq1"]):
        # SRA sample (always paired-end for now)
        accession = u["sra"]
        if sra_is_paired_end(accession):
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "sra/pe/{accession}_{group}.fastq.gz",
                        accession=accession,
                        group=["1", "2"],
                    ),
                )
            )
        else:
            return {"fq1": f"sra/se/{accession}.fastq.gz"}
    if not is_paired_end(wildcards.sample):
        return {"fq1": f"{u.fq1}"}
    else:
        return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}

def get_fastq_files(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return {"results/trimmed/{sample}_{unit}_{group}.fastq.gz".format(**wildcards)}
            
        # single end sample
        return {"results/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)}
    else:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit)]
        if pd.isna(u["fq1"]):
            # SRA sample (always paired-end for now)
            accession = u["sra"]
            if sra_is_paired_end(accession):
                if wildcards.group=='R1':
                    return {f"sra/pe/{accession}_1.fastq.gz"}
                elif wildcards.group=='R2':
                    return {f"sra/pe/{accession}_2.fastq.gz"}
            else:
                return {f"sra/se/{accession}.fastq.gz"}
        if not is_paired_end(wildcards.sample):
            return {"fq1": f"{u.fq1}"}
        else:
            if wildcards.group=='R1':
                return [f"{u.fq1}"]
            elif wildcards.group=='R2':
                return [f"{u.fq2}"]

def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

def sra_is_paired_end(sra):
    output = sp.getoutput("fastq-dump -I -X 1 -Z --split-spot "+sra+" 2>/dev/null \
        | awk '{if(NR % 2 == 1) print substr($1,length($1),1)}' \
        | uniq \
        | wc -l")
    return output=='2'

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    unit = units.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand(
            "sra/{accession}_{read}.fastq.gz", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return units.loc[wc.sample, fq].tolist()


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]
