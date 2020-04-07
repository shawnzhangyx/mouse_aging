setwd("../../analysis/snapATAC/DH/snapFiles/")
suppressPackageStartupMessages(library("SnapATAC"))

library("tictoc")
library("umap")
library("leiden")

load("DH.pool.snapATA.Frag1000.TSS10.Landmark20k.seed1.dimPC50.K20.res0.7.cluster.RData")


x.sp = x.sp[x.sp@cluster %in% c(1,2,3)]
## Run diffusion maps. 
x.sp = runDiffusionMaps(obj = x.sp,input.mat ="bmat",num.eigs = 50)
## plot dimension reduction.
plotDimReductPW(
  obj=x.sp,
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=paste(outF,params, ".PW.pdf" ,sep = ""),
  pdf.height=7,
  pdf.width=7
);

tic("runViz")
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:dims,
  method="umap",
  seed.use=10
)
toc()

# KNN
tic("runKNN")
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:dims,
  k=k
)
toc()

tic("runCluster Leiden")
x.sp=SnapATAC::runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  res = res,
  seed.use=10
)

toc()

