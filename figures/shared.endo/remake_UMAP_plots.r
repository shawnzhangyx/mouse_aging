a=read.delim("luisa_all_tissues/all_merged_40kL.cluster.meta.txt")
#a =read.delim("luisa_all_tissues/all_merged.cluster.meta.txt")
library(SnapATAC)
load("luisa_all_tissues/all_merged_40kL.cluster.RData")

b= cbind(a,x.sp@umap)
colnames(b)[c(16,17)] = c("umap.1","umap.2")
b$stage = as.character(b$stage)
b$stage = factor(b$stage,levels=c("3","10","18"))
b$cluster = factor(b$cluster)

b = b[sample(nrow(b)),]

png("umap.highlight-endo.png",height=2000,width=2000,res=300)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=cluster==27),size=0.5,alpha=0.05)+
  theme_classic(base_size=25) + #scale_color_manual(values=c("grey","
  theme(legend.position="none")
  dev.off()

  
d = b[which(b$cluster==27),]  
d2 = d[which(d$tissue_cluster %in% c("DH_13","FC_15","HT_3","HT_7","LM_6","BM_15")),]

ggplot(d2) + geom_point(aes(x=umap.1,y=umap.2,color=tissue_cluster),
  size=1,alpha=0.2) + 
  xlim(0,5) + ylim(-5,0) 

ggplot(d2) + geom_point(aes(x=umap.1,y=umap.2,color=tissue_cluster),
  size=1,alpha=0.5) +
  xlim(0,5) + ylim(-5,0) + 
  facet_wrap(.~tissue_cluster)



## 
a = read.delim("luisa_all_tissues/endo/Endo_cluster.meta.txt")
load("luisa_all_tissues/endo/Endo_.cluster.RData")

b=cbind(a,x.sp.sub@umap)
colnames(b)[c(18,19)] = c("umap.1","umap.2")

png("umap.endo-subcluster.by_tissue.png",height=2000,width=2000,res=300)
ggplot(subset(b,cluster!=5)) + geom_point(aes(x=umap.1,y=umap.2,color=tissue_cluster),size=1,alpha=0.2)+
  theme_classic(base_size=25) + #scale_color_manual(values=c("grey","
  theme(legend.position="none")
dev.off()

png("umap.endo-subcluster.facet_by_tissue.png",height=3000,width=3000,res=300)

ggplot(subset(b,cluster!=5)) + 
  geom_point(data=subset(b,cluster!=5)[,c("umap.1","umap.2")], 
                                aes(x=umap.1,y=umap.2),size=.5,alpha=0.2,color="grey") + 
  geom_point(aes(x=umap.1,y=umap.2,color=tissue_cluster),size=.5,alpha=0.2)+
  theme_classic(base_size=25) + facet_wrap(.~tissue_cluster) +
  theme(legend.position="none")
dev.off()


gglist = list()

tcs = unique(b$tissue_cluster) 

#for (tc in tcs ) {
#gglist[[tc]] =  ggplot( subset(b,cluster!=5)) + 
#  geom_point(aes(x=umap.1,y=umap.2,color=tissue_cluster==tc),size=1,alpha=0.2)+
#  theme_classic(base_size=25) + theme(legend.position="none")
#}

library(gridExtra)





  
dev.off()





