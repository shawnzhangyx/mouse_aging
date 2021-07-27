library(Seurat)
a=readRDS("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/DH_seurat_rmdoub_filtered.rds")

a$cluster_age = paste0(a$seurat_clusters,"_",a$stage)

#Idents(a) = "cluster_age"
Idents(a) = "seurat_clusters"


FeaturePlot(a,features="Ptgds",split.by="stage")


b = subset(a,seurat_clusters==1)

Idents(b) = "cluster_age"

DotPlot(b,features="Robo1")
FeaturePlot(b,features="Robo1",split.by="stage")
VlnPlot(b,features="Robo1")


FeaturePlot(a,features="Robo1",split.by="stage")
VlnPlot(subset(a,stage!=10),features="Robo1",split.by="stage")


