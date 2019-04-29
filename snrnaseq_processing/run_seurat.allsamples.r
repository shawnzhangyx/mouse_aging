library(dplyr)
library(Seurat)

# Load the dataset
#pbmc.data <- Read10X(data.dir = "JB_55/outs/filtered_feature_bc_matrix/")
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

samples = paste0("DH_",c("03_rep1","03_rep2","10_rep1","10_rep2","18_rep1","18_rep2"))

seurat.list = list()
for ( sample in samples) {
 print(sample)
 tmp.data = Read10X(data.dir = paste0("../../data/snRNA/cellranger.intron/",sample,"/outs/filtered_feature_bc_matrix/"))
 tmp = CreateSeuratObject(counts = tmp.data, project = sample, min.cells = 3, min.features = 200)
 tmp$sample= sample
 tmp$stage = sub(".._(..)_rep.","\\1",sample)
 seurat.list[[sample]] = tmp
}

data.big = merge(x=seurat.list[[1]], y = seurat.list[2:6], add.cell.ids = samples, project = "DH")
pbmc = data.big


pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
## remove unwanted cells. 
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## normalize data
#pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(object = pbmc)

pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# linear scaling of data: 
all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc)#, features = all.genes)

pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))

#pbmc <- JackStraw(object = pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
ElbowPlot(object = pbmc)


pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)

pbmc <- RunUMAP(object = pbmc, dims = 1:10)
DimPlot(object = pbmc, reduction = "umap",lable=T)

VlnPlot(object = pbmc, features = c("Mog","Apoe", "Pdgfra","C1qb",
  "Slc17a7","Gad2", "Pvalb","Sst","Npy","Vip"))

saveRDS(pbmc, file = "../../analysis/sn_rnaseq/DH_seurat.rds")

#########

#new.cluster.ids <- c("Memory CD4 T", "Naive CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
#    "NK", "DC", "Mk")
#    names(x = new.cluster.ids) <- levels(x = pbmc)
#    pbmc <- RenameIdents(object = pbmc, new.cluster.ids)
#    DimPlot(object = pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



