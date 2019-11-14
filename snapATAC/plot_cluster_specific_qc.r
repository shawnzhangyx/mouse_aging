tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/snapATAC/",tissue))
a=read.delim(paste0('snapFiles/',tissue,'.pool.snapATAC.cluster.meta.txt'))

bfiles = list.files("../../../data/snATAC/qc/bc_info_by_sample/",pattern=tissue)

ll = list()
for (idx in 1:6){
tmp=read.delim(paste0("../../../data/snATAC/qc/bc_info_by_sample/",bfiles[idx]))
tmp$sample=substr(bfiles[idx],1,10)
ll[[idx]] = tmp
}

b = do.call(rbind,ll)


b$name=paste(b$sample,b$barcodes)
a$name=paste(a$x.sp.sample,a$barcode)
b$cluster = a$x.sp.cluster[match(b$name,a$name)]

library(gridExtra)


pdf(paste0(tissue,".summary_stats.by_cluster.pdf"))
ggplot(a) + geom_boxplot(aes(x=factor(x.sp.cluster),y=TN),outlier.shape=NA) + 
  ylim(0,10000)

ggplot(a) + geom_boxplot(aes(x=factor(x.sp.cluster),y=UM),outlier.shape=NA) +
  ylim(0,10000)

ggplot(a) + geom_boxplot(aes(x=factor(x.sp.cluster),y=CM/UM),outlier.shape=NA) 

ggplot(b) + geom_boxplot(aes(x=factor(cluster),y=FIP.pc),outlier.shape=NA) 

g1 = ggplot(subset(b,!is.na(cluster))) + geom_point(aes(x=filter,y=FIP.pc,color=factor(cluster)),alpha=0.5)
g2 = ggplot(subset(b,!is.na(cluster))) + geom_histogram(aes(FIP.pc),bins=100)+coord_flip()
grid.arrange(g2,g1,widths=c(1,4))

dev.off()

