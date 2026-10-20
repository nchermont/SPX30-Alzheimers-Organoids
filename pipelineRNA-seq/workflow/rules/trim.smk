rule get_sra:
    output:
        temp("sra/pe/{accession}_1.fastq.gz"),
        temp("sra/pe/{accession}_2.fastq.gz"),
    log:
        "logs/get-sra/{accession}.log",
    params:
        extra="--skip-technical"
    threads: 6
    resources:
        load=51
    wrapper:
        "v7.0.0/bio/sra-tools/fasterq-dump"

rule get_sra_se:
    output:
        temp("sra/se/{accession}.fastq.gz"),
    log:
        "logs/get-sra/{accession}.log",
    params:
        extra="--skip-technical"
    threads: 6
    resources:
        load=51
    wrapper:
        "v7.0.0/bio/sra-tools/fasterq-dump"

rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}_{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed_cut/{sample}_{unit}_R1.fastq.gz",
        fastq2="results/trimmed_cut/{sample}_{unit}_R2.fastq.gz",
        qc="results/trimmed_cut/{sample}_{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "v7.0.0/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed_cut/{sample}_{unit}_single.fastq.gz",
        qc="results/trimmed_cut/{sample}_{unit}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-se"],
        adapters_r1=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "v7.0.0/bio/cutadapt/se"

rule trim_galore_pe:
    input:
        unpack(get_fq_trimgalore),
    output:
        fasta_fwd="results/trimmed/{sample}_{unit}_R1.fastq.gz",
        report_fwd="results/trimmed/reports/{sample}_{unit}_R1_trimming_report.txt",
        fasta_rev="results/trimmed/{sample}_{unit}_R2.fastq.gz",
        report_rev="results/trimmed/reports/{sample}_{unit}_R2_trimming_report.txt",
    threads: 1
    params:
        extra="--illumina -q 20",
    log:
        "logs/trim_galore/{sample}_{unit}.log",
    wrapper:
        "v7.0.0/bio/trim_galore/pe"

rule trim_galore_se:
    input:
        get_cutadapt_input,
    output:
        fasta="results/trimmed/{sample}_{unit}_single.fastq.gz",
        report="results/trimmed/report/{sample}_{unit}_single.fastq.gz_trimming_report.txt",
    params:
        extra="--illumina -q 20",
    log:
        "logs/trim_galore/{sample}_{unit}.log",
    wrapper:
        "v7.0.0/bio/trim_galore/se"
