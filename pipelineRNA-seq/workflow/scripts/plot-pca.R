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
svg(snakemake@output[[1]])
plotPCA(counts, intgroup = snakemake@wildcards[["variable"]])
dev.off()


svg(snakemake@output[[2]])
plotPCA(counts, intgroup = snakemake@wildcards[["variable"]]) + geom_text_repel(aes(label = name), max.overlaps = Inf)
dev.off()

