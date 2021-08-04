library(Seurat)

setwd("../../analysis/paired_tag_rna/")

a = readRDS("/projects/ren-transposon/home/chz272/transposon/13.Paired-Tag_Ageing/05.merged/02.merged_RNA/03.pileup100/2020_09_21_remove_lowqual.RDS")

meta = read.delim("../../scripts/pairedtag_rna_analysis/sample_info.txt",stringsAsFactors=F)

a$id_2 = as.numeric(substr(names(a$id),10,11))

a$sample = meta$sample[match(a$id_2,meta$ID)]
a$age = substr(a$sample,4,5)

cluster_info = read.table("../../scripts/pairedtag_rna_analysis/cluster_info.txt",header=T,stringsAsFactors=F)

a$celltype = cluster_info$celltype[match(a$seurat_clusters,cluster_info$cluster)]

cts = unique(cluster_info$celltype)

a$ct_age = factor(paste0(a$celltype,":",a$age), levels = paste0(rep(cts,each=2),":",rep(c("03","18"),length(cts))))


Idents(a) = "ct_age"


for ( ct in unique(cluster_info$celltype)){
  print(ct)
  f = FindMarkers(a, ident.1 = paste0(ct,":03"),ident.2=paste0(ct,":18"),verbose=FALSE)
  f = f[which(f$p_val <0.05),]
  write.table(f,paste0("age_diff/",ct,".03vs18.diff.txt"))
}


png("top_genes.vln.png",width=2048,height=2048)
VlnPlot(a,features=c("AC149090.1","Gm10722","Robo1","Lars2","Btaf1","Mir770","Atp2b1"),pt.size=0)
dev.off()

pdf("top_genes.dot.pdf")
DotPlot(a,features=c("AC149090.1","Gm10722","Robo1","Lars2","Btaf1","Mir770","Atp2b1"))
dev.off()

