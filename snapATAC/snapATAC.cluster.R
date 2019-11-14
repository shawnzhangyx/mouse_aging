#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("-d", "--pc_dim", default = 20, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
parser$add_argument("--tss-cutoff",type="integer", required=TRUE,help = "TSS enrichment cutoff")
parser$add_argument("--fragment-cutoff",type="integer", required=TRUE,help = "TSS enrichment cutoff")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")
library("leiden")


RDataF = args$input
dims = args$pc_dim
cpus = args$cpu
outF = args$output
TSS_cutoff = args$tss_cutoff
frag_cutoff = args$fragment_cutoff

print(TSS_cutoff)
print(frag_cutoff)
print(typeof(TSS_cutoff))
# paramters for clustering 
#dims = 20 
k=50
load(RDataF)

# apply the TSS enrichment cutoff

x.sp = x.sp[which(x.sp@metaData$TSS_enrich >= TSS_cutoff & x.sp@metaData$UQ >= frag_cutoff  )]


# new method. Create landmark cells, and run dimension reduction, and extend
# to other cells.
row.covs = log10(Matrix::rowSums(x.sp@bmat) + 1)
row.covs.dens = density(x = row.covs,bw = 'nrd', adjust = 1)
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(1)
idx.landmark.ds = sort(sample(
  x = seq(nrow(x.sp)),
  size = 10000,
  prob = sampling_prob
))
x.landmark.sp = x.sp[idx.landmark.ds, ]
x.query.sp = x.sp[-idx.landmark.ds, ]
x.landmark.sp = runDiffusionMaps(obj = x.landmark.sp,input.mat ="bmat",num.eigs = 50)
# extend from landmark cells to all cells.
x.query.sp = runDiffusionMapsExtension(obj1 = x.landmark.sp,obj2 =x.query.sp,input.mat = "bmat")
x.landmark.sp@metaData$landmark = 1
x.query.sp@metaData$landmark = 0
x.sp = snapRbind(x.landmark.sp, x.query.sp)
cat("bound\n")

## combine landmarks and query cells;
x.sp = x.sp[order(x.sp@sample),]; # IMPORTANT
rm(x.landmark.sp, x.query.sp); # free memory

plotDimReductPW(
  obj=x.sp,
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name=paste(outF, ".PW.pdf" ,sep = ""),
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

## Clustering
# R-igraph
tic("runCluster Leiden")
x.sp=SnapATAC::runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  res = 0.5,
  seed.use=10
)

toc()

plotViz(
  obj=x.sp,
  method="umap",
  main=paste(outF,"\n","TSS cutoff:",TSS_cutoff,"| Fragment cutoff:", frag_cutoff, "\n  # cells passed: ", length(x.sp@barcode)),
  point.color=x.sp@cluster,
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
  legend.add=FALSE,  pdf.file.name = paste(outF,".cluster.","TSS",TSS_cutoff,".Frag",frag_cutoff, ".pdf", sep = "")
)


## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste(outF, ".cluster.meta.txt", sep="")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF,".filterTSS", TSS_cutoff,".Frag",frag_cutoff, ".dims", dims, ".cluster.RData",sep="")
save(x.sp, file=outfname)

