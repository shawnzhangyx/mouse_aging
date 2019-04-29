library(dplyr)
library(Seurat)

m3.data = Read10X(data.dir = "../../data/snRNA/cellranger.intron/DH_03_rep1/outs/filtered_feature_bc_matrix/")
ctrl = CreateSeuratObject(counts = m3.data, project = "DH_03", min.cells = 3, min.features = 200) 
ctrl$age = "M03"
ctrl <- subset(x = ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(object = ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)

m18.data = Read10X(data.dir = "../../data/snRNA/cellranger.intron/DH_18_rep1/outs/filtered_feature_bc_matrix/")
stim = CreateSeuratObject(counts = m18.data, project = "DH_18", min.cells = 3, min.features = 200)
stim$age = "M18"
stim <- subset(x = stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(object = stim, verbose = FALSE)
stim <- FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000)


immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


DefaultAssay(object = immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(object = immune.combined, verbose = FALSE)
immune.combined <- RunPCA(object = immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

