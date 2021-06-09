setwd("../../analysis/snapATAC/")

files = list.files(pattern = "*great.red.tsv",full.names=T,recursive=T)

dat = list()

for (file in files){
  print(file)
  tissue = sub("\\.\\/(.*?)\\/.*","\\1",file)
  cluster = sub(".*great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\1",file)
  category = sub(".*great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\2",file)
  tmp=readLines(file)
  tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)

  gobp = tab[which(tab$X..Ontology=="GO Biological Process"),]
  if (length(gobp) <1) { next }
  gobp_red = gobp
  gobp_red = gobp_red[which(gobp_red$TotalGenes>5 &gobp_red$TotalGenes<1000 &
    gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1),]
  if (nrow(gobp_red) <1) { next }

#  tmp = gobp_red[,c("ID","Desc","BinomFdrQ")]
  tmp = gobp_red[,c("ID","Desc","HyperFdrQ","BinomFdrQ","ObsGenes","ExpGenes","TotalGenes")]
  dat[[file]] = data.frame(tissue=tissue,cluster=cluster,category=category,tmp)

  }

out = do.call(rbind,dat)
rownames(out)= 1:nrow(out)

meta =read.delim("../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
out$celltype = paste0(meta$Tissue,".",meta$Name)[match(paste(out$tissue,out$cluster),paste(meta$Tissue,meta$Cluster))]


up = out[which(out$category=="up" & out$BinomFdrQ < 0.001),]
#up = out[which(out$category=="up" & out$HyperFdrQ < 0.001),]

up = up[,c("celltype","Desc","BinomFdrQ")]
#up = up[,c("celltype","Desc","HyperFdrQ")]

reup = reshape(up, idvar="celltype",timevar="Desc",direction="wide")
rownames(reup) = reup$celltype
colnames(reup) = sub("BinomFdrQ.(.*?)","\\1",colnames(reup))
#colnames(reup) = sub("HyperFdrQ.(.*?)","\\1",colnames(reup))

reup$celltype=NULL
reup[is.na(reup)] = 1
reup.log = -log10(reup+0.0001)

library(pheatmap)
pdf("all_celltypes/all_celltypes.GREAT.up.pdf", width=50,height=20)
pheatmap(reup.log)
#pheatmap(reup.log,clustering_distance_cols="correlation",clustering_distance_rows="correlation")
dev.off()


down = out[which(out$category=="down" & out$BinomFdrQ < 0.001),]
#down = out[which(out$category=="down" & out$HyperFdrQ < 0.001),]

down = down[,c("celltype","Desc","BinomFdrQ")]
#down = down[,c("celltype","Desc","HyperFdrQ")]

redown = reshape(down, idvar="celltype",timevar="Desc",direction="wide")
rownames(redown) = redown$celltype
colnames(redown) = sub("BinomFdrQ.(.*?)","\\1",colnames(redown))
#colnames(redown) = sub("HyperFdrQ.(.*?)","\\1",colnames(redown))

redown$celltype=NULL
redown[is.na(redown)] = 1
redown.log = -log10(redown+0.0001)

library(pheatmap)
pdf("all_celltypes/all_celltypes.GREAT.down.pdf", width=50,height=20)
pheatmap(redown.log)
#pheatmap(redown.log,clustering_distance_cols="correlation",clustering_distance_rows="correlation")
dev.off()

