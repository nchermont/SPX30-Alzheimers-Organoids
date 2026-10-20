rule run_gsea:
    input:
        table="results/diffexp/{contrast}.diffexp.tsv",
        fpkm_path="results/fpkm/all.tsv",
    output:
        gsea_result="results/diffexp/{contrast}.gseares.RDS",
    params:
        contrast=config["diffexp"]["contrasts"],
        gsea_use_stat=config["diffexp"]["gsea_use_stat"],
    conda:
        "../envs/R_4.yaml"
    log:
        "logs/run_gsea/{contrast}.log",
    threads: 4
    resources:
        time_min=59 * 4,
        mem_mb=8192 * 7,
    script:  #
        "../scripts/run_gsea.R"