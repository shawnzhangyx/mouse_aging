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

outmetaf <- paste0("cluster1-3.meta.txt")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste("cluster1-3.subcluster.RData",sep="")
save(x.sp, file=outfname)

runMag <- function(
  obj,
  input.mat=c("gmat", "pmat", "bmat", "mmat"),
  step.size=3
){
  message("Epoch: checking the inputs ...");
  if(missing(obj)){
    stop("obj is missing");
  }else{
    if(!is(obj, "snap")){
      stop("obj is not a snap obj")
    }
  }
  ncell = nrow(obj);

  # 2. check if input matrix exists
  input.mat = match.arg(input.mat);
  if(input.mat == "bmat"){
    data.use = obj@bmat;
    peak.use = obj@feature;
  }else if(input.mat == "pmat"){
    data.use = obj@pmat;
    peak.use = obj@peak;
  }else if(input.mat == "gmat"){
    data.use = obj@gmat;
    peak.use = GenomicRanges::GRanges();
  }else{
    data.use = obj@mmat;
    peak.use = GenomicRanges::GRanges();
  }

  if((x=nrow(data.use)) == 0L){
    stop("input matrix is empty")
  }else{
    if((x=nrow(data.use)) != ncell){
      stop("input matrix has wrong number of rows with number of barcodes")
    }
  }

  cat("3! \n")
  # 3. check if the KNN graph exists
  A = obj@graph@mat;
  if((x=nrow(A)) != ncell){
    stop("affnity graph is empty, runKNN first!")
  }

  cat("4\n")
  A = A + t(A);
  A = A / Matrix::rowSums(A);
  step.size = 1;
  if(step.size > 1){
    for(i in 1:step.size){
      cat("step", i, "\n")
      cat(dim(A), "\n")
      show(head(A))
      A = A %*% A;
    }
  }
  data.use.smooth = A %*% data.use;
  cat("datause\n")
  slot(obj, input.mat) = data.use.smooth;
  return(obj)
}


