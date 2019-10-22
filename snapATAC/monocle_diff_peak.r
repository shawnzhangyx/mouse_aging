library(monocle)
library(SnapATAC)
library(chromVAR)

setwd("../../analysis/snapATAC/DH/snapFiles/")
load("DH.pool.snapATAC.cluster.RData")

x.sp = addPmatToSnap(x.sp)

x.sp.2 = makeBinary(x.sp, mat="pmat");
#data frame with barcodes in rownames and "cluster"i as single column with colname = cluster
cell_anno_df=data.frame(cluster=paste("cluster", x.sp@cluster, sep = ""), row.names=paste(x.sp@sample, x.sp@barcode, sep = "_"))
#binary cell x peak matrix - if using peak x cell - matrix remove transpose on line 16
ncounts=x.sp.2@pmat
#remove blacklist regions
#black_list_1=read.table("/mnt/tscc/rraviram/hg38/hg38.blacklist.bed", stringsAsFactors = FALSE)
#human=read.table("/mnt/tscc/rraviram/scATAC-Seq/Human_brain_high_regions.bed", stringsAsFactors = FALSE)
#black_list=rbind(black_list_1, human[,1:3])
#black_list.gr = GRanges(
# black_list[,1],
# IRanges(black_list[,2], black_list[,3])
#);
#idy_peak = queryHits(findOverlaps(x.sp@peak, black_list.gr));
#ncounts=t(ncounts[,-idy_peak])

#column names is cell barcodes
ncounts = t(ncounts)
colnames(ncounts) = rownames(cell_anno_df)
#peaks data frame
ygi_df=data.frame(x.sp@peak)[,1:3]
rownames(ygi_df)=paste(ygi_df[,1],"_", ygi_df[,2],"_", ygi_df[,3], sep = "")
#rownames is peak naems
rownames(ncounts)=rownames(ygi_df)

pd <- new("AnnotatedDataFrame",data=cell_anno_df)
fd <- new("AnnotatedDataFrame",data=ygi_df)
cds <- newCellDataSet(ncounts,phenoData = pd, featureData = fd,expressionFamily=binomialff(),
                     lowerDetectionLimit=1)
pData(cds)$Size_Factor = 1
pData(cds)$Age = substr(rownames(pData(cds)),4,5)
pData(cds)$Rep = substr(rownames(pData(cds)),7,10)

saveRDS(cds,"DH.pool.snapATAC.monocle_CDS.RData")
#diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~cluster", cores=5)


### subsetting cluster 6
c6 = cds[,rownames(pData(cds))[which(pData(cds)$cluster=="cluster6")]]
### differential binding test. 
diff_test_res <- differentialGeneTest(c6,
                    fullModelFormulaStr = "~Age",cores=20)

diff_test_res = diff_test_res[order(diff_test_res$pval),]

