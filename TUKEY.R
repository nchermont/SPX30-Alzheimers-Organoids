library(readxl)
# Table format: first column (label) with conditions. Following columns: accession  
# Replace the file path, changing '\' to '/'  
tabela <- read_excel("/Users/namoreira/Desktop/Book12.xlsx")

pca <- prcomp(tabela[,-1], scale.=TRUE)
summary(pca)

library(ggfortify)
autoplot(pca, data=tabela, colour='label')

testes <- data.frame()
for(i in 2:ncol(tabela)){
  a <- data.frame(label = tabela[,1], value = tabela[,i])
  protname <- colnames(a)[2]
  colnames(a) <- c("label","value")
  model <- aov(value~label, data=a)
  tuk <- TukeyHSD(model, conf.level=.95)
  ab <- as.data.frame(tuk$label)
  ab$comparacao <- rownames(ab)
  rownames(ab) <- NULL
  ab$proteina <- protname
  testes <- rbind(testes,ab)
}

library(tidyverse)
testes <- testes %>%
  relocate(proteina,comparacao)
print(testes)
library(openxlsx)
#choose path
setwd("/Users/namoreira/Downloads")
write.xlsx(testes, file="Tukey_RESULTS_genotype.xlsx")