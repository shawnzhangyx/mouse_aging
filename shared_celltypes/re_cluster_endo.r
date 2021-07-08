#Endo:
# DH.13
# FC.15
# LM.6
# HT.7
# HT.3
# BM.15



library(SnapATAC)
#load("/projects/ps-renlab/lamaral/projects/Aging/all_tissues/all_merged_40kL.cluster.RData")

load("/mnt/tscc/lamaral/projects/Aging/all_tissues/endo/Endo_.cluster.RData")


tab = cbind(x.sp.sub@metaData,x.sp.sub@umap,x.sp.sub@cluster)

colnames(tab)[c(17,18,19)] = c("umap1","umap2","subcluster")

ggplot(tab) + geom_point(aes(x=umap1,y=umap2,color=tissue_cluster),alpha=0.5)+
  facet_wrap(.~tissue_cluster)


ggplot(tab) + geom_point(aes(x=umap1,y=umap2,color=tissue_cluster),alpha=0.5)

ggplot(tab) + geom_point(aes(x=umap1,y=umap2,color=cluster),alpha=0.5) +
  facet_wrap(.~cluster)

ggplot(tab) + geom_point(aes(x=umap1,y=umap2,color=subcluster),alpha=0.5)

ggplot(tab) + geom_point(aes(x=umap1,y=umap2,color=subcluster),alpha=0.5) + 
  facet_wrap(.~subcluster)

  table(tab$tissue_cluster,tab$subcluster)


######### calculate the doublet scores. 
tab$is5 = tab$subcluster==5
#st = tab[which(tab$subcluster==5),]

db = read.delim("/mnt/tscc/lamaral/projects/Aging/HT/doublets/all_doub_meta_bmat.txt")

tab$db = db$doublet_scores[match( paste0(tab$sample,tab$barcode), paste0(db$sample, db$barcode))]

tab2 = tab[which(tab$tissue=="HT"),]

ggplot(tab2) + geom_histogram(aes(db)) + facet_grid(is5~.)

ggplot(tab2) + geom_point(aes(x=umap1,y=umap2,color=db))

ggplot(tab2) + geom_histogram(aes(db)) + facet_grid(tissue_cluster~.)

