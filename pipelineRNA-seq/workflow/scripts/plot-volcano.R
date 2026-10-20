log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("KernSmooth")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

contrast_config <- snakemake@config[["diffexp"]][["contrasts"]][[
    snakemake@wildcards[["contrast"]]
]]

# basic case of contrast specification, see:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
if (length(contrast_config) == 2 && typeof(contrast_config) == "list") {
  if (
    # check for existence contrast's variable_of_interest to
    # provide useful error message
    !(contrast_config[["variable_of_interest"]] %in%
      names(snakemake@config[["diffexp"]][["variables_of_interest"]])
    )
  ) {
      cli_abort(
        c(
                "config.yaml: All variable_of_interest entries under `diffexp: contrasts:`",
          " " = "must also exist under `diffexp: variables_of_interest:`.",
          "x" = "Could not find variable_of_interest: {contrast_config[['variable_of_interest']]}",
          " " = "It was not among the `diffexp: variables_of_interest:`",
          " " = "{names(snakemake@config[['diffexp']][['variables_of_interest']])}",
          "i" = "Are there any typos in the contrasts' `variable_of_interest:` entries?"
        )
      )
  }
  contrast <- c(
    contrast_config[["variable_of_interest"]],
    contrast_config[["level_of_interest"]],
    snakemake@config[["diffexp"]][["variables_of_interest"]][[
      contrast_config[["variable_of_interest"]]
    ]][["base_level"]]
  )
# more complex contrast specification via list(c(), c()), see ?results docs of
# the DESeq2 package and this tutorial (plus the linked seqanswers thread):
# https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
} else if (
    length(contrast_config) == 1 &&
    typeof(contrast_config) == "character"
  ) {
  contrast <- d <- eval(parse(text = contrast_config))
}

res <- results(dds, contrast = contrast)

alpha <- 0.01 # Threshold on the p-value

# par(mfrow=c(1,2))

# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res$sig <- -log10(res$padj)
sum(is.infinite(res$sig))

res[is.infinite(res$sig),"sig"] <- 350
# View(resultDESeq2[is.na(resultDESeq2$pvalue),])

# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(res$pvalue)
# sum(genes.to.plot)
range(res[genes.to.plot, "log2FoldChange"])

## Volcano plot of adjusted p-values
cols <- densCols(res$log2FoldChange, res$sig)
cols[res$pvalue ==0] <- "purple"
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

svg(snakemake@output[[1]])

plot(res$log2FoldChange, 
     res$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=res$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
#gn.selected <- abs(res$log2FoldChange) > 2 & res$padj < alpha 
#if(length(res$log2FoldChange[gn.selected])>0){
#     text(res$log2FoldChange[gn.selected],
#     -log10(res$padj)[gn.selected],
#     lab=rownames(res)[gn.selected ], cex=0.6)
#}


dev.off()
