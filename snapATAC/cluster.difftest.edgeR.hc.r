setwd("../../analysis/snapATAC/DH/age_diff_edgeR.snap")
a=read.csv("4.edger.txt")


diff = a[which(a$PValue<0.01),c(2:7)]

cor.mat = cor(diff)
hc.sample = hclust(as.dist(1-cor.mat))
cor.mat = cor(t(diff))
hc.peak = hclust(as.dist(1-cor.mat))

diff = diff[hc.peak$order,hc.sample$order]
rownames(diff) = 1:nrow(diff)

melted = melt(as.matrix(diff))
ggplot(melted) +geom_tile(aes(x=Var2,y=Var1,fill=log10(value+0.1))) 




diff.scale = t(apply(diff,1,scale))
colnames(diff.scale) = colnames(diff)
melted = melt(as.matrix(diff.scale))

ggplot(melted) +geom_tile(aes(x=Var2,y=Var1,fill=value)) +
  scale_fill_gradient2(high="red",low="blue",mid="white")
