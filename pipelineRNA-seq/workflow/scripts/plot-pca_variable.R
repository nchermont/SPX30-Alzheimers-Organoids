log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library(ggplot2)
library(ggrepel)

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])



# obtain normalized counts
counts <- rlog(dds, blind=FALSE)

counts.sub <- counts[ , counts[[snakemake@wildcards[["variable"]]]] %in% c(snakemake@wildcards[["value"]]) ]

svg(snakemake@output[[1]])
if(ncol(counts.sub)>2)
    plotPCA(counts.sub, intgroup = snakemake@wildcards[["variable"]])
dev.off()


svg(snakemake@output[[2]])
if(ncol(counts.sub)>2)
    plotPCA(counts.sub, intgroup = snakemake@wildcards[["variable"]]) + geom_text_repel(aes(label = name), max.overlaps = Inf)
dev.off()

