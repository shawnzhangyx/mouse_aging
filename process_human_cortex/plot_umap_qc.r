#tissue="DH"
tissue=commandArgs(trailing=T)[1]
library(viridis)

setwd(paste0("../../analysis/human_cortex/",tissue))

a=read.delim(paste0(tissue,".pool.barcode.meta_info.txt"))

pdf("cluster.UMAP.QC.pdf",height=5,width=8)
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

## plot the TSS percent on the UMAP
#ggplot(a) + geom_point(aes(umap.1,umap.2,color=TSS.pc),size=.1,alpha=0.5) +
#  scale_color_viridis(direction=1) +
#  scale_color_gradientn(colors=c("red","blue","darkblue"),values=c(0,0.3,1)) +
#  facet_grid(rep~stage) + theme_bw()




ggplot(a) + geom_boxplot(aes(sample,TSS_enrich))

ggplot(a) + geom_boxplot(aes(factor(cluster),TSS_enrich))
#ggplot(a) + geom_violin(aes(factor(cluster),TSS.pc))


ggplot(a)  + geom_boxplot(aes(factor(cluster),TSS_enrich)) +
  facet_grid(~sample)

#ggplot(a) + geom_histogram(aes(TSS_enrich),bins=300)

#ggplot(a) + geom_histogram(aes(TSS_enrich),bins=100) + 
#  facet_grid(rep~stage) +
#  geom_vline(xintercept=7,linetype='dashed') 



dev.off()
