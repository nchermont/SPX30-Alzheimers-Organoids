rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}_{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}_{unit}.log",
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
    threads: 30
    wrapper:
        "v7.0.0/bio/star/align"

rule samtools_index:
    input:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools_index/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
    threads: 30  # This value - 1 will be sent to -@
    wrapper:
        "v7.0.0/bio/samtools/index"