library(Seurat)
library(doParallel)
library(ggrepel)
library(plotly)
library(gridExtra)
library(edgeR)

tissue=commandArgs(trailing=T)[1]
#tissue="DH"

setwd("../../analysis/repeat_analysis/")
## read the matrix
pbmc = readRDS(paste0(tissue,"_seurat.rds"))

## read the class and family inforamtion.
class= read.delim("../../analysis/repeat_analysis/repeat_class_family.txt",header=F)
class$V1 = gsub("_","-",class$V1)
class$V1 = gsub("\\(","-",class$V1)
class$V1 = gsub("\\)","-",class$V1)

## extract the counts from the object. 
mat = as.matrix(pbmc[["RNA"]]@counts)


# mat = sweep(mat,2,colSums(mat),'/') # normalize the counts per cell. 

#mat = mat[which(rownames(mat) %in% c("LSU-rRNA-Hsa","GSAT-MM","L1MdA-III")),]
sample = pbmc$sample
cluster = pbmc$cluster
sb = names(pbmc$sample)
#meta = read.delim("../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.meta.txt")
#meta$sb = paste(meta$sample,meta$barcode,sep="_")

# calculate the sum of counts by sample. 
mat_sample = t(rowsum(t(mat),group=sample))
mat_sample = mat_sample[order(-mat_sample[,1]),]

mat1 = data.frame(sweep(mat_sample,2,colSums(mat_sample),'/'))
mat1$name = rownames(mat1)
mat1 = mat1[-which(mat1$name=="Others"),]

#g1 = ggplot(mat1, aes(x=DH_03_rep1,y=DH_03_rep2,text=name)) + geom_point() +
#  geom_abline(intercept=0,slope=1,color='red')#+scale_x_log10() +scale_y_log10()
# ggplotly(g1)

#g2 = ggplot(mat1, aes(x=DH_10_rep1,y=DH_10_rep2,text=name)) + geom_point() +
#  geom_abline(intercept=0,slope=1,color='red')  #+scale_x_log10() +scale_y_log10()
#ggplotly(g2)

y= DGEList(counts=mat_sample[-which(rownames(mat_sample) %in% c("Others","GSAT-MM")),])
#y= DGEList(counts=mat_sample)
y = calcNormFactors(y)

avg = as.data.frame(cpm(y))/1e4
avg$name = rownames(avg)
avg$class = class$V2[match(avg$name,class$V1)]

x = 1; y=2
gglist = list()
for (x in 1:6){
  for ( y in 1:6){
    gglist[[length(gglist)+1]] = ggplot(avg) + geom_point(aes_string(x=colnames(avg)[x],y=colnames(avg)[y]),size=0.1)+
      geom_abline(intercept=0,slope=1,color='red') + 
      xlab(colnames(avg)[x]) + ylab(colnames(avg)[y]) +
      theme_bw()
  }
}

pdf(paste0(tissue,".repeat.all_sample.scatter.pdf"),height=10,width=10)
grid.arrange(grobs=gglist,nrow=6)
dev.off()


#ggplotly(ggplot(avg) + geom_point(aes(x=DH_10_rep1,y=DH_10_rep2,label=name)))

#pdf(paste0(tissue,".repeat.10_rep1.vs.rep2.pdf"),height=10,width=10)
#ggplot(avg) + geom_point(aes(x=DH_10_rep1,y=DH_10_rep2,color=class=="LINE/L1")) + 
#  geom_text_repel(data=avg[which(grepl("L1Md",avg$name) & avg$DH_10_rep2/avg$DH_10_rep1>1.5),],
#                  aes(x=DH_10_rep1,y=DH_10_rep2,label=name)) +
#  geom_text_repel(data=subset(avg,grepl("rRNA",avg$name)),aes(x=DH_10_rep1,y=DH_10_rep2,label=name) )+
#  geom_abline(intercept=0,slope=1,color='red') 
#dev.off() 


#pdf(paste0(tissue,".repeat.03_rep1.vs.rep2.pdf"),height=10,width=10)
#ggplot(avg) + geom_point(aes(x=DH_03_rep1,y=DH_03_rep2,color=class=="LINE/L1")) + 
#  geom_text_repel(data=subset(avg,grepl("rRNA",avg$name)),aes(x=DH_03_rep1,y=DH_03_rep2,label=name) )+
#  geom_abline(intercept=0,slope=1,color='red') 
#dev.off() 

#ggplot(avg) + geom_point(aes(x=DH_10_rep1,y=DH_10_rep2,color=class=="SINE/Alu")) + 
#  geom_text_repel(data=avg[which(grepl("L1Md",avg$name) & avg$DH_10_rep2/avg$DH_10_rep1>1.5),],
#                  aes(x=DH_10_rep1,y=DH_10_rep2,label=name)) +
#  geom_text_repel(data=subset(avg,grepl("rRNA",avg$name)),aes(x=DH_10_rep1,y=DH_10_rep2,label=name) )+
#  geom_abline(intercept=0,slope=1,color='red') 

## calculate the sum of counts by sample and cluster. 
#mat_cluster_sample = t(rowsum(t(mat),group=paste0(cluster,".",sample)))
#mat_cluster_sample = mat_cluster_sample[order(-mat_cluster_sample[,1]),]

#mat1 = data.frame(sweep(mat_cluster_sample,2,colSums(mat_cluster_sample),'/'))
#mat1$name = rownames(mat1)
#mat1 = mat1[-which(mat1$name=="Others"),]

#g1 = ggplot(mat1, aes(x=X2.DH_03_rep1,y=X11.DH_03_rep1,text=name)) + geom_point() +
#  geom_abline(intercept=0,slope=1,color='red') + scale_x_log10() +scale_y_log10()

#ggplotly(g1)

