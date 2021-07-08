library(SnapATAC)
#load("/projects/ps-renlab/lamaral/projects/Aging/all_tissues/all_merged_40kL.cluster.RData")

load("/mnt/tscc/lamaral/projects/Aging/all_tissues/endo/Endo_.cluster.RData")


tab = cbind(x.sp.sub@metaData,x.sp.sub@umap,x.sp.sub@cluster)

colnames(tab)[c(17,18,19)] = c("umap1","umap2","subcluster")

a=read.delim("../HT/HT.pool.barcode.meta_info.txt")


a$db = paste(a$sample, a$barcode) %in% paste(tab$sample,tab$barcode)[which(tab$subcluster==5)]

a = a[order(a$db),]
ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=db))

ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=db)) + 
  facet_wrap(.~cluster) 

