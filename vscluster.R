install.packages("BiocManager")
BiocManager::install("S4Arrays")
library(S4Arrays)
library(vsclust)
library(openxlsx)
#paste the file path changing the bar to "/" 
dat <- read.xlsx("/Users/namoreira/Library/CloudStorage/GoogleDrive-nchmoreira0@gmail.com/My Drive/POSDOC UCSD/Experimentos/space project/spx30 data/proteomics AD /IDEP data_log10.xlsx", rowNames = TRUE)
#condition names, add and exclude as needed, change 5 for the number of replicates in each condition
nomes_col <- c(paste("Ground_WTAD", 1:3, sep=""), 
               paste("GROUND_APP", 1:3, sep=""), 
               paste("GROUND_PSEN1", 1:3, sep=""),
               paste("GROUND_PSEN2", 1:2, sep=""),
               paste("ISS_WTAD", 1:3, sep=""),
               paste("ISS_APP", 1:3, sep=""),
               paste("ISS_PSEN1", 1:3, sep=""),
               paste("ISS_PSEN2", 1:3, sep=""))
colnames(dat) <- nomes_col
dat <- dat[,order(as.numeric(substr(colnames(dat), nchar(colnames(dat)), nchar(colnames(dat)))))]
dat[] <- lapply(dat, function(x) as.numeric(as.character(x)))

library(vsclust)
#### running statistical analysis and estimation of individual variances
# dat, NUMBER OF REPLICATES, NUMBER OF CONDITIONS, FALSE, TRUE
statOut <- PrepareForVSClust(dat, 3, 8, FALSE, TRUE)

dat <- statOut$dat
Sds <- dat[,ncol(dat)]
cat(paste("Features:",nrow(dat),"\nMissing values:",
          sum(is.na(dat)),"\nMedian standard deviations:",
          round(median(Sds,na.rm=TRUE),digits=3)))
## Write output into file 
write.csv(statOut$statFileOut,
          paste("","Prot","vsclustResult_SPX31AD.csv",sep=""))


ClustInd <- estimClustNum(dat, maxClust=11, scaling="standardize", cores=4)
PreSetNumClustVSClust <- 0
PreSetNumClustStand <- 0
#### Use estimate cluster number or use own
if (PreSetNumClustVSClust == 0){
  PreSetNumClustVSClust <- optimalClustNum(ClustInd)
}
if (PreSetNumClustStand == 0){
  PreSetNumClustStand <- optimalClustNum(ClustInd, method="FCM")
}
#### Visualize
estimClust.plot(ClustInd)

library(matrixStats)
#### Run clustering (VSClust and standard fcm clustering
ClustOut <- runClustWrapper(dat, 
                            PreSetNumClustVSClust, 
                            NULL, 
                            VSClust=TRUE, 
                            scaling="standardize",
                            cores=1)

Bestcl <- ClustOut$Bestcl
VSClust_cl <- Bestcl
## Write clustering results (VSClust)
write.csv(data.frame(cluster=Bestcl$cluster,
                     ClustOut$outFileClust,
                     isClusterMember=rowMaxs(Bestcl$membership)>0.5,
                     maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), 
          paste("Prot", 
                "FCMVarMResults_VSclust_SPX30AD", 
                Sys.Date(), 
                ".csv", 
                sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, 
          paste("Prot",
                "FCMVarMResultsCentroids_SPX30AD",
                Sys.Date(), 
                ".csv", 
                sep=""))

ClustOut <- runClustWrapper(dat, PreSetNumClustStand, NULL, VSClust=FALSE, 
                            scaling="standardize", cores=1)
Bestcl <- ClustOut$Bestcl
## Write clustering results (standard fcm)
write.csv(data.frame(cluster=Bestcl$cluster,
                     ClustOut$outFileClust,
                     isClusterMember=rowMaxs(Bestcl$membership)>0.5,
                     maxMembership=rowMaxs(Bestcl$membership),
                     Bestcl$membership), 
          paste("Prot", 
                "FCMResults_CLASSIC_SPX30AD", 
                Sys.Date(), 
                ".csv", 
                sep=""))
## Write coordinates of cluster centroids
write.csv(Bestcl$centers, paste("Prot",
                                "FCMResultsCentroids_CLASSIC_SPX30AD", 
                                Sys.Date(),
                                ".csv", 
                                sep=""))