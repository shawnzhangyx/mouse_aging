setwd("../../analysis/snapATAC/all_celltypes/")

a = data.frame(fread("all_celltypes.counts"))
b = read.delim("../../../aging_share/figures/celltype_annotation.txt",sep="\t")

cnts = a[,-(1:6)]
tissue = sub("(..).*metacell_(.*).sorted.bam","\\1",colnames(cnts))
cluster = sub("(..).*metacell_(.*).sorted.bam","\\2",colnames(cnts))

names = paste0(b$Tissue,".",b$Cluster,".",b$Name)[match(paste(tissue,cluster),paste(b$Tissue,b$Cluster))]
colnames(cnts) = names
# remove low quality 
cnts = cnts[,-which(paste(tissue,cluster) %in% paste(b$Tissue,b$Cluster)[is.na(b$Clade)])]

library(edgeR)
y = DGEList(counts=cnts)
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)

#mds = plotMDS(y)
#dat = data.frame(name=colnames(y$counts),x=mds$x,y=mds$y,tissue=substr(colnames(y$counts),1,2))
#library(ggrepel)
#ggplot(dat) + geom_label_repel(aes(x=x,y=y,label=name,color=tissue))


cpms = cpm(y)
std = apply(cpms,1,sd)
ave = apply(cpms,1,mean)
std.norm = std/ave
cpms2 = cpms[which(std>=quantile(std,0.50)),]
#cpms2 = cpms[which(std.norm>=quantile(std.norm,0.80)),]

log.cpm = log2(cpms2+0.01)

dist.sp = cor(cpms2,method="spearman")
hc.sp = hclust(d=as.dist(1-dist.sp))
dist.log = cor(log.cpm)
hc.log = hclust(d=as.dist(1-dist.log))


pdf("all_celltypes.hierarchical_clustering.pdf",width=10,height=15)
par(mar=c(2,10,2,10))
plot(as.dendrogram(hc.sp),horiz=T,main="Spearman")
plot(as.dendrogram(hc.log),horiz=T,main="log.Cpm")

dev.off()

#library(pheatmap)
#pheatmap(dist)


