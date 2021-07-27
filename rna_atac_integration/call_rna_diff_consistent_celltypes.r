library(Seurat)
library(ggplot2)

setwd("../../analysis/rna_atac_integration")

pbmc = readRDS("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/DH_seurat_rmdoub_filtered.rds")
meta =read.table("../../scripts/rna_atac_integration/rna_cell_type.consistent.txt",header=T)
#pbmc$seurat_clusters.stage = paste0(pbmc$seurat_clusters,"_",pbmc$stage)
#Idents(object = pbmc) <- "seurat_clusters.stage"
pbmc$celltype = meta$celltype[match(pbmc$seurat_clusters,meta$cluster)]

pbmc$celltype.stage = paste0(pbmc$celltype,"_",pbmc$stage)

celltypes = unique(pbmc$celltype)
celltypes = celltypes[!is.na(celltypes)]
library(doParallel)
registerDoParallel(cores=20)

foreach(i=celltypes) %dopar% {
#for(i in celltypes) {
  print(i)
  Idents(object = pbmc) <- "celltype.stage"

  f = FindMarkers(object = pbmc, ident.1 = paste(i,"_03", sep = ""), ident.2 = paste(i,"_18",sep = ""),min.pct = .001,logfc.threshold = .01)
#  f = FindMarkers(object = pbmc, ident.1 = paste(i,"_03", sep = ""), ident.2 = paste(i,"_18",sep = ""),test.use="LR")
#,verbose = FALSE)

  head(x = f, n = 15)
  write.table(f[which(f$p_val<1.1),], paste("RNA_diff/",i,"_03vs18_nocutoff.txt", sep = ""), sep = "\t", quote = F)
}

