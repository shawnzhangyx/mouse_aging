
setwd("../../analysis/snapATAC/")
#",tissue,"/age_diff_edgeR.snap/motif.homer"))

#files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

files = list.files(pattern = "*great.red.tsv",full.names=T,recursive=T)


dat = list()
for (file in files){
  print(file)
  
  tissue = sub("\\.\\/(.*?)\\/.*","\\1",file)
  cluster = sub("\\.\\/..\\/age_diff_edgeR.snap\\/great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\1",file)
  category = sub("\\.\\/..\\/age_diff_edgeR.snap\\/great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\2",file)
  tmp=readLines(file)
  tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)

  gobp = tab[which(tab$X..Ontology=="GO Biological Process"),]
  if (length(gobp) <1) { next }
  gobp_red = gobp
  gobp_red = gobp_red[which(gobp_red$TotalGenes>5  & gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1),]
  tmp = gobp_red[1:50,c("ID","Desc","BinomFdrQ","HyperFdrQ","ObsGenes","TotalGenes","RegionFoldEnrich")]

 # tmp = read.delim(file)
  dat[[file]] = data.frame(tissue=tissue,cluster=cluster,category=category,tmp)
  }

out = do.call(rbind,dat)
rownames(out)= 1:nrow(out)

meta =read.delim("../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
out$celltype = paste0(meta$Tissue,".",meta$Cluster,".",meta$Name)[match(paste(out$tissue,out$cluster),paste(meta$Tissue,meta$Cluster))]

outp = out[,c(1,2,11,3,4,5,6,10)]
rownames(outp) = NULL
write.csv(outp, "all_celltypes/all_celltypes.great.gobp.csv",row.names=F)


up = out[which(out$category=="up" & out$BinomFdrQ<0.05 & 
          out$HyperFdrQ<0.01 & out$RegionFoldEnrich>=2),]

#up$celltype = paste(up$tissue,up$cluster)
up = up[,c("celltype","Desc","BinomFdrQ")]
reup = reshape(up, idvar="celltype",timevar="Desc",direction="wide")
rownames(reup) = reup$celltype
colnames(reup) = gsub(" ",".",sub("BinomFdrQ.(.*?)","\\1",colnames(reup)))
reup$celltype=NULL
reup[is.na(reup)] = 1
reup.log = -log10(reup+0.0001)
reup.log = t(reup.log)

library(pheatmap)
pdf("all_celltypes/all_celltypes.Great.up.pdf", width=20,height=60)
#pheatmap(reup.log)
pheatmap(reup.log)
#,clustering_distance_cols="correlation",clustering_distance_rows="correlation") 
dev.off()



