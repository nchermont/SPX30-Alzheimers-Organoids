log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("pheatmap")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

df <- as.data.frame(colData(dds)[c(snakemake@wildcards[["variable"]])])
ntd <- normTransform(dds)

#plot top 20 genes
top20 <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
rownames(df) <- colnames(assay(ntd)[top20,])
svg(snakemake@output[["top20"]])
pheatmap(assay(ntd)[top20,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

#plot top 100 genes
top100 <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
rownames(df) <- colnames(assay(ntd)[top100,])
svg(snakemake@output[["top100"]])
pheatmap(assay(ntd)[top100,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

res <- results(dds)

#plot all deg genes
alldeg_genes <- rownames(res)[!is.na(res$padj) & res$padj < 0.05]
rownames(df) <- colnames(assay(ntd)[alldeg_genes,])
svg(snakemake@output[["all"]])
if(nrow(assay(ntd)[alldeg_genes,])>=2){
    pheatmap(assay(ntd)[alldeg_genes,], cluster_rows=TRUE, show_rownames=FALSE,
            cluster_cols=TRUE, annotation_col=df)
}
dev.off()

#plot up regulated genes
upregulated_genes <- rownames(res)[!is.na(res$padj) & res$log2FoldChange > 0 & res$padj < 0.05]
rownames(df) <- colnames(assay(ntd)[upregulated_genes,])
svg(snakemake@output[["up"]])
if(nrow(assay(ntd)[upregulated_genes,])>=2){
    pheatmap(assay(ntd)[upregulated_genes,], cluster_rows=TRUE, show_rownames=FALSE,
            cluster_cols=TRUE, annotation_col=df)
}
dev.off()

#plot down regulated genes
downregulated_genes <- rownames(res)[!is.na(res$padj) & res$log2FoldChange < 0 & res$padj < 0.05]
rownames(df) <- colnames(assay(ntd)[downregulated_genes,])
svg(snakemake@output[["down"]])
if(nrow(assay(ntd)[downregulated_genes,])>=2){
    pheatmap(assay(ntd)[downregulated_genes,], cluster_rows=TRUE, show_rownames=FALSE,
            cluster_cols=TRUE, annotation_col=df)
}
dev.off()