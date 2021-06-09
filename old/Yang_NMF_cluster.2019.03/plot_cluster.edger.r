tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]
setwd(paste0("../../analysis/Yang_NMF_method/",tissue,"/R",rank))
system("mkdir cluster_feature")
a=read.delim(paste0(tissue,".counts"),skip=1)
b =read.delim(paste0(tissue,".counts.summary"))
cnts = a[,-c(1:6)]
total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])

grps = paste0("C",sub(".*metacell_(.{1,2}\\....rep.).bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)


library(edgeR)

y = DGEList(cnts)
y$sample$lib.size=total
y = calcNormFactors(y)

coord = plotMDS(y)

dat = data.frame(name=grps,x=coord$x,y=coord$y)
dat$cluster = sub("(C.*)\\.(..)\\.(rep.)","\\1",dat$name)
dat$age = sub("(C.*)\\.(..)\\.(rep.)","\\2",dat$name)
dat$rep = sub("(C.*)\\.(..)\\.(rep.)","\\3",dat$name)

library(ggrepel)
pdf("MDS_plot.pdf",width=10,height=10)
ggplot(dat) +  geom_point(aes(x=x,y=y,color=cluster)) + 
  geom_label_repel(aes(x=x,y=y,color=cluster,label=name))
ggplot(dat) + geom_point(aes(x=x,y=y,color=age)) +
  geom_label_repel(aes(x=x,y=y,color=age,label=name))
ggplot(dat) + geom_point(aes(x=x,y=y,color=rep)) +
  geom_label_repel(aes(x=x,y=y,color=rep,label=name))
dev.off()

