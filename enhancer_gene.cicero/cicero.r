library(cicero)
library(SnapATAC)
library(Matrix)

tissue =commandArgs(trailing=T)[1]
snapFile=commandArgs(trailing=T)[2]

setwd(paste0("../../analysis/cicero_results/",tissue))
load(snapFile)

barcode = paste0(x.after.sp@sample,".",x.after.sp@barcode)
peaks = x.after.sp@peak
mat = x.after.sp@pmat

# read and transpose the matrix
indata = t(x.after.sp@pmat)
indata@x[indata@x>0]<-1

# format cell info
cellinfo = data.frame(cells=paste0(x.after.sp@sample,".",x.after.sp@barcode))
row.names(cellinfo) <- cellinfo$cells

# format peak info
peakinfo = data.frame(x.after.sp@peak@seqnames,x.after.sp@peak@ranges)[,1:3]
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
## make input CDS. 
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=VGAM::binomialff(),
                            lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)
#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

# 
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)

#dimension reduction. 
umap_coords<- x.after.sp@umap
umap_coords<-as.matrix(umap_coords)
rownames(umap_coords)<-row.names(pData(input_cds))
cicero_cds<-make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

saveRDS(cicero_cds, file=paste0(tissue,".cicero.rds"))
##
library(doParallel)
registerDoParallel(cores=4)

mm10.genome = read.table("../../../annotations/mm10.chrom.sizes",stringsAsFactors=F)

output= foreach(chr=mm10.genome$V1[1:21]) %dopar% {
sample_genome <- subset(mm10.genome, V1 == chr )
conns <- run_cicero(cicero_cds, sample_genome) # Takes a few minutes to run
conns
}
conns = do.call(rbind,output)
#head(conns)

write.csv(conns,paste0(tissue,".cicero_conns.csv"))



