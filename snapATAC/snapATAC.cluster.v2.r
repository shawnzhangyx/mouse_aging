#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-i", "--input", required=TRUE, help="load preprocessed RData")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# parameters for barcode selection.
parser$add_argument("--tss-cutoff",type="integer",default=7, help = "TSS enrichment cutoff")
parser$add_argument("--fragment-cutoff",type="integer",default=500,help = "Fragment cutoff")
# parameters for dimension reduction. 
parser$add_argument("-d", "--pc_dim", default = 50, help="num of PCA dims used for clustering [default %(default)s]")
parser$add_argument("--seed", default=1, type="integer",help="Random seed for landmark selection.")
parser$add_argument("--landmark_num",default=10000, type="integer", help="Number of landmark cells.") 
# paramters for clustering. 
parser$add_argument("--knn_k", default=20, type="integer",help="K for KNN.")
parser$add_argument("--leiden_res", default=0.7, type="double",help="resolution parameter for leiden.")
# parameters for plot. 
parser$add_argument("--marker_genes", type="character",help="path to marker gene list")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library("tictoc")
library("umap")
library("leiden")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(GenomicRanges))



RDataF = args$input
cpus = args$cpu
outF = args$output
TSS_cutoff = args$tss_cutoff
frag_cutoff = args$fragment_cutoff
landmark_num = args$landmark_num
dims = args$pc_dim
seed = args$seed
k = args$knn_k
res = args$leiden_res
marker_genes = args$marker_genes

params = paste0(".Frag",frag_cutoff,".TSS",TSS_cutoff,".Landmark",as.integer(landmark_num/1000),"k.seed",seed,".dimPC",dims,".K",k,".res",res)


load(RDataF)

# apply the TSS enrichment cutoff
#x.sp = x.sp[which(x.sp@sample == "DH_03_rep1")]
x.sp = x.sp[which(x.sp@metaData$TSS_enrich >= TSS_cutoff & x.sp@metaData$UQ >= frag_cutoff  )]

print(paste0("Number of Cells:",nrow(x.sp@metaData)))

# new method. Create landmark cells, and run dimension reduction, and extend
# to other cells.
row.covs = log10(Matrix::rowSums(x.sp@bmat) + 1)
row.covs.dens = density(x = row.covs,bw = 'nrd', adjust = 1)
sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
set.seed(seed)
idx.landmark.ds = sort(sample(
  x = seq(nrow(x.sp)),
  size = landmark_num, #default 10,000
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

## Clustering
# R-igraph
tic("runCluster Leiden")
x.sp=SnapATAC::runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  res = res,  
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
  legend.add=FALSE,  pdf.file.name = paste(outF,params,".cluster.UMAP.pdf", sep = "")
)


## Create chromatin lanscape and identify cis-elements for each cluster seperately.

outmetaf <- paste0(outF,params,".meta.txt")
outmetamx <- cbind(x.sp@sample, x.sp@metaData, x.sp@cluster, x.sp@umap)
write.table(outmetamx, outmetaf, row.names=F, col.names=T, sep="\t", quote=F)

outfname = paste(outF,params, ".cluster.RData",sep="")
save(x.sp, file=outfname)

genes = read.table("gencode.vM16.gene.bed")
genes.gr = GRanges(genes[,1],
    IRanges(genes[,2], genes[,3]), name=genes[,4]
     );

marker.genes = readLines(marker_genes)
# c("Mog","Apoe","Pdgfra","C1qb","Slc17a7","Gad2","Nov","Prox1","Foxd1")

genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];

x.sp = createGmatFromMat(
    obj=x.sp,
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=10
  );
#
x.sp = scaleCountMatrix(
    obj=x.sp,
    cov=x.sp@metaData$UQ + 1,
    mat="gmat",
    method = "RPM"
  );


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


  # 4. smooth
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


x.sp = runMag(
  obj=x.sp,
    input.mat="gmat",
      step.size=1
      )

genes = genes.sel.gr$name
pdf(paste(outF,params, ".markers_UMAP.pdf", sep = ""), height = 7, width = 8)
par(mfrow = c(3,3),oma = c(0, 0, 2, 0))
for(i in genes){
  plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@gmat[, i],
    method="umap",
    main=paste(i),
    point.size=0.1,
    point.shape=19,
    down.sample=10000,
    quantiles=c(.01,.99),
    pdf.file.name = NULL
  )}
mtext(outF, outer = TRUE, cex = 1.5)
dev.off()

##### plot UMAP QC. 
a = read.delim(outmetaf)

a$stage = substr(a$sample,4,5)
a$rep = substr(a$sample,7,10)


tab = table(a$x.sp.cluster,a$sample)
tab = sweep(tab,2,colSums(tab),'/')

melted = melt(tab)


pdf(paste0(outF,params,".cluster.UMAP.QC.pdf"),height=5,width=8)

## plot the number of reads on the UMAP
ggplot(a) + geom_point(aes(umap.1,umap.2,color=log10(UM)),size=0.1,alpha=0.5) +
  scale_color_viridis(direction=1) +
#  scale_color_gradientn(colors=c("red","blue","darkblue"),values=c(0,0.5,1)) #+
  facet_grid(rep~stage) + theme_bw()


## plot the TSS enrichment core on the UMAP.
ggplot(a) + geom_point(aes(umap.1,umap.2,color=TSS_enrich),size=.1,alpha=0.5) +
#  scale_color_viridis(direction=1) +
  scale_color_gradientn(colors=c("red","blue","darkblue"),values=c(0,0.1,1)) +
  facet_grid(rep~stage) + theme_bw()

ggplot(a) + geom_boxplot(aes(sample,TSS_enrich))
ggplot(a) + geom_boxplot(aes(factor(x.sp.cluster),TSS_enrich))
ggplot(a)  + geom_boxplot(aes(factor(x.sp.cluster),TSS_enrich)) +
  facet_grid(~sample)
## plot cell type fraction. 
ggplot(melted) +
  geom_col(aes(factor(Var1),value,fill=factor(substr(Var2,4,5))
    ,color=factor(substr(Var2,7,10))),position="dodge") +
  xlab("Cell Types") + ylab("Fraction") +
  guides(fill=guide_legend(title="Age Group"))  +
  theme_bw()

dev.off()



