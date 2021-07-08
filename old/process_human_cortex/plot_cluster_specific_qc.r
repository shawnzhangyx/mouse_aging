tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/human_cortex/",tissue))

a=read.delim(paste0(tissue,".pool.barcode.meta_info.txt"))


library(gridExtra)


pdf(paste0(tissue,".summary_stats.by_cluster.pdf"))
#ggplot(a) + geom_boxplot(aes(x=factor(cluster),y=TN),outlier.shape=NA) + 
#  ylim(0,10000)

ggplot(a) + geom_boxplot(aes(x=factor(cluster),y=UM),outlier.shape=NA) +
  ylim(0,10000)

ggplot(a) + geom_boxplot(aes(x=factor(cluster),y=CM/UM),outlier.shape=NA) 

dev.off()


