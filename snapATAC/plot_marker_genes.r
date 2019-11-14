library(GenomicRanges);
library(SnapATAC);


outF = "DH"
setwd("../../analysis/snapATAC/DH/snapFiles/")

load("DH.pool.snapATAC.filterTSS10.Frag500.dims20.cluster.RData")

genes = read.table("../../../../scripts/snapATAC/gencode.vM16.gene.bed")
genes.gr = GRanges(genes[,1], 
    IRanges(genes[,2], genes[,3]), name=genes[,4]
     );

marker.genes = c("Mog","Apoe","Pdgfra","C1qb","Slc17a7","Gad2","Nov","Pink8", "Lpl", "Prox1","Foxd1")

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

## the following script from Luisa:
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
pdf(paste(outF, "_markers_UMAP.pdf", sep = ""), height = 7, width = 8)
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

