library(dplyr)
library(Seurat)

pbmc.data <- Read10X(data.dir = "../../analysis/repeat_analysis/DH/DH_03_rep1")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "DH", min.cells = 1, min.features = 1)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = pbmc), 10)


plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

d=read.delim("../../analysis/snapATAC/DH/DH.pooled.barcode.cluster.stage.rep.txt")
d = d[which(d$stage==3 & d$rep=="rep1"),]

pbmc$cluster = d$cluster[match(rownames(pbmc@meta.data),d$barcode)]
Idents(object = pbmc) <- "cluster"

VlnPlot(object = pbmc, features=c("MuRRS4-int","MURVY-int"))


f = FindMarkers(object = pbmc, ident.1 = 4,
    verbose = FALSE)

VlnPlot(object = pbmc, features=rownames(f)[1:5])

