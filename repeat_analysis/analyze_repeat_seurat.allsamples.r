library(dplyr)
library(Seurat)
library(stringr)

tissue="DH"
tissue=commandArgs(trailing=T)[1]

samples = paste0(tissue,"_",c("03_rep1","03_rep2","10_rep1","10_rep2","18_rep1","18_rep2"))
seurat.list = list()
for ( sample in samples) {
 print(sample)
 tmp.data <- Read10X(data.dir = paste0("../../analysis/repeat_analysis/",tissue,"/",sample))
 tmp <- CreateSeuratObject(counts = tmp.data, project = sample, min.cells = 1, min.features = 1)
  tmp$sample= sample
  tmp$stage = sub(".._(..)_rep.","\\1",sample)
  seurat.list[[sample]] = tmp
  }

pbmc = merge(x=seurat.list[[1]], y = seurat.list[2:6], add.cell.ids = samples, project = "DH")

#pbmc[["percent.simple"]] <- PercentageFeatureSet(object = pbmc, pattern = "n$")
#summary(pbmc$percent.simple)

pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = pbmc), 10)


plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))

d=read.delim("../../analysis/snapATAC/DH/DH.pooled.barcode.cluster.stage.rep.txt")
d$names = paste0(tissue,"_",str_pad(d$stage,2,'left',"0"),"_",d$rep,"_",d$barcode)

pbmc$cluster = d$cluster[match(rownames(pbmc@meta.data),d$names)]
pbmc = subset( pbmc, cluster>0)


### Make clusters. 
#pbmc <- ScaleData(object = pbmc)#, features = all.genes)
#pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
#ElbowPlot(object = pbmc)
#pbmc <- FindNeighbors(object = pbmc, dims = 1:30)
#pbmc <- FindClusters(object = pbmc, resolution = 0.1)
#pbmc <- RunUMAP(object = pbmc, dims = 1:30)
#DimPlot(object = pbmc, reduction = "umap",label=T)


saveRDS(pbmc, file = paste0("../../analysis/repeat_analysis/",tissue,"_seurat.rds"))

