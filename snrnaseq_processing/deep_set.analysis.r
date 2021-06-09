library(Seurat)

#a =readRDS("/projects/ps-renlab/yanxiao/projects/mouse_aging/aging_share/luisa_DH_rna/RNA_analysis_after_doublet_removal/DH_filtered.RDS")
a = readRDS("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/DH_seurat_rmdoub_filtered.rds")

VlnPlot(a,features=c("Cdkn2a"))

FeaturePlot(a, features = c("Cdkn2a","Slc17a7"),split.by="stage")
FeaturePlot(a, features = c("Cdkn2a","AC149090.1"),split.by="stage")


VlnPlot(subset(a,stage %in% c("18","10")), features="AC149090.1",split.by="stage",pt.size=0)
VlnPlot(a, features="AC149090.1",split.by="stage")

