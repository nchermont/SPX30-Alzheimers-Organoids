library(dplyr)
library(tidyverse)
library(tcltk)
data <- read.xlsx("/Users/namoreira/Library/CloudStorage/GoogleDrive-nchmoreira0@gmail.com/My Drive/POSDOC UCSD/Experimentos/space project/spx30 data/proteomics AD /ANOVA TABLE.xlsx")
nrep <- ncol(data)
#factors should be in columns 2 and 3, ex: column 2 cell type, column 3 sex
factor1 <- names(data)[2]
factor2 <- names(data)[3]
pvals <- data.frame ("p_values" = c(factor1, factor2 , paste(factor1, ":", factor2, "interaction")))
i=4

print(ncol(data))
print(names(data))

for(i in 4:nrep) {
  data[,i] <- as.numeric(data[,i])
  # extract the protein accession ID from the dataframe to use it in ANOVA.
  prot_acc <- names(data)[i]
  # Calculates 2-way ANOVA: creates the formula in the variable "form" then uses in in aov
  form <- as.formula(paste(prot_acc, " ~ ", factor1," * ", factor2))
  res.aov2 <- aov(form, data = data)
  # get just the p-values from the result (first element, and from that element, the fifth column)
  pvalues.res.aov2 <- summary(res.aov2)[[1]][5]
  # delete row "Residuals" (row4), but keep dataframe structure and row names
  pvalues.res.aov2 <- pvalues.res.aov2[-c(4), , drop=FALSE]
  # rename the column to match the protein ID
  names(pvalues.res.aov2)[1] <- prot_acc
  # Add a column to pvalues.res.aov2 identical to the p_values column in pvals for matching during merge
  pvalues.res.aov2 <- mutate(pvalues.res.aov2, "p_values" = c("Genotype", "Environment" , "Genotype:Environment"))
  #merges the calculated p-values into the dataframe that will contain all p-values for all proteins
  pvals <- merge(x=pvalues.res.aov2,
                 y=pvals,
                 by="p_values",
                 all=TRUE)
}
print(pvals)
#choose path
setwd("/Users/namoreira/Library/CloudStorage/GoogleDrive-nchmoreira0@gmail.com/My Drive/POSDOC UCSD/Experimentos/space project/spx30 data/proteomics AD")
write.xlsx(pvals, file ="Result_Anova2way_spx30_protein.xlsx")