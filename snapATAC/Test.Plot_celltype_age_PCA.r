setwd("../../analysis/snapATAC/DH")

a=read.delim("DH.peaks.counts",skip=1)

library(edgeR)

cnts = a[,-c(1:6)]
names = sub(".*metacell_(.*).sorted.bam","C\\1",colnames(cnts))
colnames(cnts) = names

y = DGEList(counts=cnts)
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)

## MDS Plot
mds = plotMDS(y)
cluster = sub("(C.*)\\...\\.rep.","\\1",names)
age = sub("C.*\\.(..)\\.rep.","\\1",names)
rep = sub("C.*\\...\\.(rep.)","\\1",names)

dat = data.frame(names=names,cluster=cluster,age=age,rep=rep,x=mds$x,y=mds$y)

ggplot(dat) + geom_text(aes(x=x,y=y,label=cluster,color=cluster)) 
ggplot(dat) + geom_text(aes(x=x,y=y,label=cluster,color=age))


ggplot(dat) + geom_text(aes(x=x,y=y,label=rep,color=age)) +  
  facet_wrap(~cluster) +
  scale_color_brewer(palette="RdBu")

ggplot(subset(dat, age!="03" | rep!="rep1")) + geom_text(aes(x=x,y=y,label=rep,color=age)) + 
  facet_wrap(~cluster)


cpms = cpm(y)


pca =  prcomp(t(cpms),center=T)
plot(pca$x,type='n')
text(pca$x,labels=meta$HPV.ID)


dat = data.frame(names=names,cluster=cluster,age=age,rep=rep,x=pca$x[,1],y=pca$x[,2])

