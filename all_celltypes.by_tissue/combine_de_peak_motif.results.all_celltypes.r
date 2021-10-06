
setwd("../../analysis/snapATAC/")
#",tissue,"/age_diff_edgeR.snap/motif.homer"))

files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

files = files[grep("motif.homer.bg",files)]

dat = list()
for (file in files){
  tissue = sub("\\.\\/(.*?)\\/.*","\\1",file)
  cluster = sub(".*\\/(.*)\\.(.*)\\.homer.*","\\1",file)
  category = sub(".*\\/(.*)\\.(.*)\\.homer.*","\\2",file)
  tmp = read.delim(file)[1:50,]
  dat[[file]] = data.frame(tissue=tissue,cluster=cluster,category=category,tmp[,c(1,2,3,5)])
  }

out = do.call(rbind,dat)
rownames(out)= 1:nrow(out)

meta =read.delim("../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
out$celltype = paste0(meta$Tissue,".",meta$Name)[match(paste(out$tissue,out$cluster),paste(meta$Tissue,meta$Cluster))]
out = out[-which(out$celltype %in% paste0(meta$Tissue,".",meta$Name)[is.na(meta$Clade)]),]

outp = out[,c(1,2,8,3,4,5,6,7)]
outp = outp[which(outp$P.value<1),]

write.csv(outp, "all_celltypes/all_celltypes.motifs.csv",row.names=F)
write.csv(outp[which(outp$q.value..Benjamini.<0.05),], "all_celltypes/all_celltypes.motifs.q0.05.csv",row.names=F)

up = outp[which(outp$category=="up" & outp$q.value..Benjamini.<0.05),]

#up$celltype = paste(up$tissue,up$cluster)
up = up[,c("celltype","Motif.Name","q.value..Benjamini.")]
reup = reshape(up, idvar="celltype",timevar="Motif.Name",direction="wide")
rownames(reup) = reup$celltype
colnames(reup) = sub("q.value..Benjamini..(.*?)\\/.*","\\1",colnames(reup))
reup$celltype=NULL
reup[is.na(reup)] = 1
#reup = reup[,which(colSums(reup<0.05) >1)]

reup.log = -log10(reup+1e-5)

library(pheatmap)
pdf("all_celltypes/all_celltypes.Motif.up.bg.pdf", width=35,height=15)
#pheatmap(reup.log)
pheatmap(reup.log,
#  clustering_distance_cols="correlation", 
#clustering_method = "single"
#  clustering_distance_rows="correlation"
)
dev.off()

down = outp[which(outp$category=="down" & outp$q.value..Benjamini.<0.05),]

#down$celltype = paste(down$tissue,down$cluster)
down = down[,c("celltype","Motif.Name","q.value..Benjamini.")]
redown = reshape(down, idvar="celltype",timevar="Motif.Name",direction="wide")
rownames(redown) = redown$celltype
colnames(redown) = sub("q.value..Benjamini..(.*?)\\/.*","\\1",colnames(redown))
redown$celltype=NULL
redown[is.na(redown)] = 1
#redown = redown[,which(colSums(redown<0.05) >1)]

redown.log = -log10(redown+1e-5)

library(pheatmap)
pdf("all_celltypes/all_celltypes.Motif.down.bg.pdf", width=35,height=15)
#pheatmap(redown.log)
pheatmap(redown.log,
  #clustering_distance_cols="correlation",
  #clustering_distance_rows="correlation"
  )
dev.off()

# up 
up.cnt = table(up$Motif.Name)
up.agg = aggregate(celltype~Motif.Name,up,function(vec){paste(vec,collapse=",")})
up.agg$count = up.cnt[match(up.agg$Motif.Name, names(up.cnt))]
up.agg = up.agg[order(-up.agg$count),]

# down 
down.cnt = table(down$Motif.Name)
down.agg = aggregate(celltype~Motif.Name,down,function(vec){paste(vec,collapse=",")})
down.agg$count = down.cnt[match(down.agg$Motif.Name, names(down.cnt))]
down.agg = down.agg[order(-down.agg$count),]

write.csv(up.agg,"frequency_up_Motifs.csv")
write.csv(down.agg,"frequency_down_Motifs.csv")

