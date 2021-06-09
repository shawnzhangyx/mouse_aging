
setwd("../../analysis/snapATAC/")
#",tissue,"/age_diff_edgeR.snap/motif.homer"))

files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

files = files[grep("motif.homer.bg",files)]

dat = list()
for (file in files){
  tissue = sub("\\.\\/(.*?)\\/.*","\\1",file)
  cluster = sub(".*\\/(.*)\\.(.*)\\.homer.*","\\1",file)
  category = sub(".*\\/(.*)\\.(.*)\\.homer.*","\\2",file)
  tmp = read.delim(file)
  dat[[file]] = data.frame(tissue=tissue,cluster=cluster,category=category,tmp[,c(1,2,3,5)])
  }

out = do.call(rbind,dat)
rownames(out)= 1:nrow(out)

meta =read.delim("../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
out$celltype = paste0(meta$Tissue,".",meta$Name)[match(paste(out$tissue,out$cluster),paste(meta$Tissue,meta$Cluster))]

up = out[which(out$category=="up" & out$q.value..Benjamini.<0.1),]

#up$celltype = paste(up$tissue,up$cluster)
up = up[,c("celltype","Motif.Name","q.value..Benjamini.")]
reup = reshape(up, idvar="celltype",timevar="Motif.Name",direction="wide")
rownames(reup) = reup$celltype
colnames(reup) = sub("q.value..Benjamini..(.*?)\\/.*","\\1",colnames(reup))
reup$celltype=NULL
reup[is.na(reup)] = 1
reup.log = -log10(reup+0.0001)

library(pheatmap)
pdf("all_celltypes/all_celltypes.Motif.up.pdf", width=20,height=10)
#pheatmap(reup.log)
pheatmap(reup.log,clustering_distance_cols="correlation",clustering_distance_rows="correlation")
dev.off()

down = out[which(out$category=="down" & out$q.value..Benjamini.<0.1),]

#down$celltype = paste(down$tissue,down$cluster)
down = down[,c("celltype","Motif.Name","q.value..Benjamini.")]
redown = reshape(down, idvar="celltype",timevar="Motif.Name",direction="wide")
rownames(redown) = redown$celltype
colnames(redown) = sub("q.value..Benjamini..(.*?)\\/.*","\\1",colnames(redown))
redown$celltype=NULL
redown[is.na(redown)] = 1
redown.log = -log10(redown+0.0001)

library(pheatmap)
pdf("all_celltypes/all_celltypes.Motif.down.pdf", width=20,height=10)
#pheatmap(redown.log)
pheatmap(redown.log,clustering_distance_cols="correlation",clustering_distance_rows="correlation")
dev.off()


