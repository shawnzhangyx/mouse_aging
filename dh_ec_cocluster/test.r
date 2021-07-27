library(SnapATAC)

load("/mnt/tscc/lamaral/projects/Aging/DH_EC/EC_DH.cluster.meta.RData")


a = data.frame(x.sp@metaData,x.sp@umap,x.sp@cluster)


png("alltissue.png",height=2000,width=2000)
ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),alpha=0.2)
dev.off()

png("alltissue.facet.png",height=2000,width=2000)
ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),alpha=0.2) +
  facet_wrap(.~tissue)
dev.off()

png("DH.label.png",height=2000,width=2000)
ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=transfer_DH_ann)) +
  facet_wrap(.~transfer_DH_ann)
dev.off()

ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=transfer_DH_ann==7),alpha=0.2)


b = subset(a,umap.1<0 & umap.1> -10 & umap.2 > -12 & umap.2 < -2)

ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=transfer_DH_ann==7),alpha=0.2)

ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=transfer_DH_ann==7),alpha=0.2) +
  facet_wrap(.~transfer_brain_ann)

ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=tissue),alpha=0.2)  


ggplot(b) + geom_point(aes(x=umap.1,y=umap.2,color=transfer_brain_ann),alpha=0.2) + 
  facet_wrap(.~tissue)

png("DH.7.black.ITHGL.red.png",height=2000,width=2000)
ggplot() +  geom_point(data=a, aes(x=umap.1,y=umap.2),color="grey",alpha=0.1) +  
  geom_point(data=subset(a,transfer_DH_ann==7),aes(x=umap.1,y=umap.2),alpha=0.2) + 
  geom_point(data=subset(a,transfer_brain_ann=="ITHGL"),aes(x=umap.1,y=umap.2),color="red",alpha=0.03)
dev.off()

png("DH.7.black.CA1GL.red.png",height=2000,width=2000)
ggplot() +  geom_point(data=a, aes(x=umap.1,y=umap.2),color="grey",alpha=0.1) +
  geom_point(data=subset(a,transfer_DH_ann==7),aes(x=umap.1,y=umap.2),alpha=0.2) +
  geom_point(data=subset(a,transfer_brain_ann=="CA1GL"),aes(x=umap.1,y=umap.2),color="red",alpha=0.03)
dev.off()

png("DH.7.black.NPGL.red.png",height=2000,width=2000)
ggplot() +  geom_point(data=a, aes(x=umap.1,y=umap.2),color="grey",alpha=0.1) +
  geom_point(data=subset(a,transfer_DH_ann==7),aes(x=umap.1,y=umap.2),alpha=0.2) +
  geom_point(data=subset(a,transfer_brain_ann=="NPGL"),aes(x=umap.1,y=umap.2),color="red",alpha=0.03)
dev.off()


#NPGL ITHGL CA1GL

table(b$x.sp.cluster,b$transfer_DH_ann)
# 8 599
# 14 480
# 17 178
# 20 690
table(b$x.sp.cluster,b$transfer_brain_ann)[c(8,14,17,20),]

# 8 L6bGL(5924) NPGL(1573)
# 14 AmyExN(530) CA1GL(8996) 
# 17 AmyExN(1109) Ent_neur(5342) ITHGL(2195)
# 20 ITHGL(5857)

