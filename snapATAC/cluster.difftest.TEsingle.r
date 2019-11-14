setwd("../../analysis/snapATAC/DH/snapFiles/")
load("DH.pool.snapATAC.cluster.RData")

pmat = x.sp@pmat
sample = x.sp@sample
cluster = x.sp@cluster
idx = which( sample %in% c("DH_03_rep2",
    "DH_18_rep2") & cluster == 6)

grps = paste(cluster,sample,sep=":")[idx]
cnts = Matrix::t(pmat[idx,])
rownames(cnts) = x.sp@peak$name

cnts  = cnts[which(rowSums(cnts)>median(rowSums(cnts))),]


library(DEsingle)
library(BiocParallel)

# Set the parameters and register the back-end to be used
param <- MulticoreParam(workers = 10, progressbar = TRUE)
register(param)

# Detecting the DE genes in parallelization with 18 cores
results <- DEsingle(counts = cnts, group = factor(grps), parallel = TRUE, BPPARAM = param)

peakS = rowSums(cnts)
results$peakSum = peakS[match(rownames(results),rownames(cnts))]

write.csv(results, "../DESingle.03-18.rep2.csv")

