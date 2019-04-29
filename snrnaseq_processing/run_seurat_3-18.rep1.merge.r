library(dplyr)
library(Seurat)

m3.data  = Read10X(data.dir = "../../data/snRNA/cellranger.intron/DH_03_rep1/outs/filtered_feature_bc_matrix/")
m3 = CreateSeuratObject(counts = m3.data, project = "DH_03", min.cells = 3, min.features = 200)
m3$age = "M03"

m18.data  = Read10X(data.dir = "../../data/snRNA/cellranger.intron/DH_18_rep1/outs/filtered_feature_bc_matrix/")
m18 = CreateSeuratObject(counts = m18.data, project = "DH_18", min.cells = 3, min.features = 200)
m18$age = "M18"

data.big = merge(x=m3, y = m18, add.cell.ids = c("M03","M18"), project = "DH")
pbmc = data.big

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

pbmc <- ScaleData(object = pbmc)#, features = all.genes)

pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)

pbmc <- RunUMAP(object = pbmc, dims = 1:10)
DimPlot(object = pbmc, reduction = "umap",split.by="age",label=T)

VlnPlot(object = pbmc, features = c("Mog","Apoe", "Pdgfra","C1qb",
  "Slc17a7","Gad2", "Pvalb","Sst","Npy","Vip"))


