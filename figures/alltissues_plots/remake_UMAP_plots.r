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

png("umap.cluster.png",height=3000,width=3000,res=300)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=cluster),size=0.5,alpha=0.05)+
  theme_classic() +
  theme(legend.position="none")
  dev.off()

png("umap.cluster.facet.png",height=3000,width=3000,res=300)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=cluster),size=0.5,alpha=0.05)+
  theme_classic() +
  facet_wrap(.~cluster) +
  theme(legend.position="none")
  dev.off()



#pdf("umap.tissue.pdf",height=8,width=8)
png("umap.tissue.png",height=720,width=720)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),size=0.5,alpha=0.05)+
  theme_classic() +
  theme(legend.position="none")
dev.off()

png("umap.tissue.facet.png",height=720,width=720)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),size=0.5,alpha=0.05)+
  facet_wrap(.~tissue) + 
  theme_classic()
  dev.off()


#b =b[order(b$stage),]
png("umap.age.png",height=1500,width=1500,res=300)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=stage),size=0.1,alpha=0.05)+
  scale_color_brewer(palette ="YlGnBu")+
  theme_classic() +
  theme(legend.position="none")

  dev.off()

png("umap.age.facet.png",height=720,width=1440)
ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=stage),size=1,alpha=0.05)+
  scale_color_brewer(palette ="YlGnBu")+ 
  facet_wrap(.~stage) +
  theme_classic()
    dev.off()

##endothelial cells.
b=b[which(b$cluster==27),]

ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),size=3,alpha=0.2)+
  theme_classic() + scale_x_continuous(limits=c(0,5)) + 
  scale_y_continuous(limits=c(-5,0))#+ 
#  theme(legend.position="none")

ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=stage),size=3,alpha=0.2)+
  theme_classic() + scale_x_continuous(limits=c(0,5)) +
  scale_y_continuous(limits=c(-5,0)) + 
  facet_wrap(.~stage) 

  


