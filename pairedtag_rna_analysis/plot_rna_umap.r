library(Seurat)

setwd("../../analysis/paired_tag_rna/")

a = readRDS("/projects/ren-transposon/home/chz272/transposon/13.Paired-Tag_Ageing/05.merged/02.merged_RNA/03.pileup100/2020_09_21_remove_lowqual.RDS")

meta = read.delim("../../scripts/pairedtag_rna_analysis/sample_info.txt",stringsAsFactors=F)

a$id_2 = as.numeric(substr(names(a$id),10,11))

a$sample = meta$sample[match(a$id_2,meta$ID)]
a$age = substr(a$sample,4,5)

cluster_info = read.table("../../scripts/pairedtag_rna_analysis/cluster_info.txt",header=T,stringsAsFactors=F)

a$celltype = cluster_info$celltype[match(a$seurat_clusters,cluster_info$cluster)]

Idents(a) = "celltype"

png("umap.rna.png",width=720,height=560)
DimPlot(a)
dev.off()



a$celltype.age = paste0(a$celltype,"_",a$age)
a$rep = sub(".*_(Rep.)","\\1",a$sample)
a$celltype.age.rep = paste0(a$celltype.age,"_",a$rep)

Idents(a) = "celltype.age.rep"

b = subset(a, celltype== "OGC")

dot = DotPlot(b,features="Apoe")


