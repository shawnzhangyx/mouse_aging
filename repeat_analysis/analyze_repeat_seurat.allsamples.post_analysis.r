library(Seurat)
library(doParallel)
library(ggrepel)
tissue="DH"
tissue=commandArgs(trailing=T)[1]

setwd("../../analysis/repeat_analysis/")
## read the matrix
pbmc = readRDS(paste0(tissue,"_seurat.rds"))

## read the class and family inforamtion.
class= read.delim("../../analysis/repeat_analysis/repeat_class_family.txt",header=F)
class$V1 = gsub("_","-",class$V1)
class$V1 = gsub("\\(","-",class$V1)
class$V1 = gsub("\\)","-",class$V1)


#Idents(object = pbmc) <- "cluster"
#f = FindMarkers(object = pbmc, ident.1 = 4,
#    verbose = FALSE)
#VlnPlot(object = pbmc, features=rownames(f)[1:5])

Idents(object = pbmc) <- "stage"
avg <- log1p(x = AverageExpression(pbmc, verbose = FALSE)$RNA)
avg$family = class$V3[match(rownames(avg),class$V1)]
avg$class = class$V2[match(rownames(avg),class$V1)]
avg$name=rownames(avg)
avg = subset(avg,name!="Others")

pdf(paste0(tissue,".repeat_03.vs.18.all.pdf"),height=10,width=10)
ggplot(avg) + geom_point(aes(x=`03`,y=`18`,color=family)) +
  geom_label_repel(data=subset(avg,family=="SINE"),aes(x=`03`,y=`18`,label=name))+
  geom_abline(intercept=0,slope=1) + 
  ggtitle("Combined")
#ggplot(avg) + geom_point(aes(x=`03`,y=`18`,color=family=="SINE")) +
#  ggtitle("Combined")
ggplot(avg) + geom_point(aes(x=`03`,y=`18`,color=class)) +
  ggtitle("Combined")
dev.off()

#registerDoParallel(cores=max(pbmc$cluster)/2)

avg.list=list()
Idents(object = pbmc) <- "cluster"
for (i in 1:max(pbmc$cluster)){
print(i)
c5 = subset(pbmc, cluster==i)
Idents(c5) <- "stage"
avg.c5 <- log1p(x = AverageExpression(c5, verbose = FALSE)$RNA)
avg.c5$class = class$V2[match(rownames(avg.c5),class$V1)]
avg.c5$family = class$V3[match(rownames(avg.c5),class$V1)]
avg.c5 = subset(avg.c5,rownames(avg.c5)!="Others")

print(head(avg.c5))
avg.list[[i]] = avg.c5
}

pdf(paste0(tissue,".repeat_03.vs.18.clusters.pdf"),height=8,width=8)
for (cluster in 1:max(pbmc$cluster)){
print(cluster)
print(ggplot(avg.list[[cluster]]) + geom_point(aes(x=`03`,y=`18`,color=class=="LINE/L1")) +   geom_abline(intercept=0,slope=1) +
  ggtitle(paste0("cluster-",cluster)) )
  }
dev.off()


pdf(paste0(tissue,".repeat_03.vs.18.clusters.family.pdf"),height=8,width=8)
for (cluster in 1:max(pbmc$cluster)){
print(cluster)
print(ggplot(avg.list[[cluster]]) + geom_point(aes(x=`03`,y=`18`,color=family)) +
  geom_abline(intercept=0,slope=1) +
  ggtitle(paste0("cluster-",cluster)) )
  }
dev.off()

