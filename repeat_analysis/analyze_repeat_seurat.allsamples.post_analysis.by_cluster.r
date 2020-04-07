library(Seurat)
library(doParallel)
library(ggrepel)
library(plotly)
library(gridExtra)
library(edgeR)

tissue=commandArgs(trailing=T)[1]
tissue="FC"

setwd("../../analysis/repeat_analysis/")
## read the matrix
pbmc = readRDS(paste0(tissue,"_seurat.rds"))


## extract the counts from the object. 
mat = as.matrix(pbmc[["RNA"]]@counts)

#mat = mat[which(rownames(mat) %in% c("LSU-rRNA-Hsa","GSAT-MM","L1MdA-III")),]
sample = pbmc$sample
cluster = pbmc$cluster
sb = names(pbmc$sample)


## calculate the sum of counts by sample and cluster. 
mat_cluster_sample = t(rowsum(t(mat),group=paste0(cluster,".",sample)))
mat_cluster_sample = mat_cluster_sample[order(-mat_cluster_sample[,1]),]

cpms = sweep(mat_cluster_sample,2,colSums(mat_cluster_sample),"/")*1e6

for (clu in sort(unique(cluster)) ) {

idx = which(colnames(mat_cluster_sample) %in% 
  paste0(clu, ".", tissue, c("_03_rep1","_03_rep2",
  "_10_rep1","_10_rep2",   "_18_rep1","_18_rep2")))

y = DGEList(counts=mat_cluster_sample[,idx])
y = calcNormFactors(y)
grps = c("03","03","10","10","18","18")
design = model.matrix(~0+grps)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast =c(-1,0,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]
out = out[-which(rownames(out) %in% c("Others")),]
print(c(clu, length(which(out$fdr<0.05))))
}


ggplotly(ggplot(out) + geom_point(aes(x=`3.DH_03_rep2`, y=`3.DH_18_rep2`,label=rownames(out))) + 
          geom_abline(intercept=0,slope=1)
         
         )


mat1 = data.frame(sweep(mat_cluster_sample,2,colSums(mat_cluster_sample),'/'))
mat1$name = rownames(mat1)
mat1 = mat1[-which(mat1$name=="Others"),]

g1 = ggplot(mat1, aes(x=X2.DH_03_rep1,y=X11.DH_03_rep1,text=name)) + geom_point() +
  geom_abline(intercept=0,slope=1,color='red') + scale_x_log10() +scale_y_log10()

ggplotly(g1)




