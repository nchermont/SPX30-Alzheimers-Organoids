# Tables and figures of the paper, see docs/paper_figures.md


rule paper_figures:
    input:
        counts="results/counts/all.tsv",
        samples=config["samples"],
        annotation=config["paper"]["annotation"],
    output:
        directory("results/paper_figures"),
    log:
        "logs/paper_figures.log",
    conda:
        "../envs/paper_figures.yaml"
    script:
        "../scripts/paper_figures.R"
