library(harmony)
library(SnapATAC)
library(leiden)

dims = 20
k=20
res=0.7
RDataF = commandArgs(trailing=T)[1]
  outF = commandArgs(trailing=T)[2]
params=".harmony"

load(RDataF)

x.after.sp = runHarmony(
  obj=x.sp, 
  eigs.dim=1:dims, 
  meta_data=x.sp@sample # sample index
  );

x.after.sp = runViz(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:dims,
  method="umap",
  seed.use=10
)


x.after.sp = runKNN(
  obj=x.after.sp,
  eigs.dims=1:dims,
  k=k
)

x.after.sp=SnapATAC::runCluster(
  obj=x.after.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  res = res,
  seed.use=10
)


plotViz(
  obj=x.after.sp,
  method="umap",
  main=paste(outF,"\n  # cells passed: ", length(x.after.sp@barcode)),
  point.color=x.after.sp@cluster,
  point.size=.5,
  point.shape=19,
  point.alpha=0.8,
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE,  pdf.file.name = paste(outF,params,".cluster.UMAP.pdf", sep = "")
)

outfname = paste(outF,params, ".cluster.RData",sep="")
save(x.after.sp, file=outfname)

outmetaf <- paste0(outF,params,".meta.txt")
outmetamx <- cbind(x.after.sp@sample, x.after.sp@metaData, x.after.sp@cluster, x.after.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

