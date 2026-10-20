log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("cli")
library(stringr)
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

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
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]

norm_counts = counts(dds, normalized=T)

res <- results(
  dds,
  contrast = contrast,
  parallel = parallel
)


# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial, see:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
res <- lfcShrink(
  dds,
  contrast = contrast,
  res = res,
  type = "ashr"
)

# sort by p-value
res <- res[order(res$padj), ]


res <- data.frame(
    "gene" = rownames(res),
    res
  )

head(res)
conditions <- names(snakemake@config[["diffexp"]][["variables_of_interest"]])
for (condition in conditions)
{
  for (i in unique(coldata[[condition]]))
  {
      mean <- data.frame(mean=rowMeans(norm_counts[,c(coldata[coldata[[condition]]==i,]$sample_name)]))
      mean <- cbind(rownames(mean), mean)
      rownames(mean) <- NULL
      colnames(mean) <- c("gene", paste("baseMean_", condition, "_", i, sep=""))
      head(mean)
      res<-merge(res, mean)
  }
}
res<-merge(res, norm_counts,by.x="gene", by.y='row.names')


# TODO explore IHW usage


write.table(
  res,
  file = snakemake@output[["table"]],
  row.names = FALSE,
  sep = "\t"
)
