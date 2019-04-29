library(Seurat)


setwd("../../analysis/sn_rnaseq")

pbmc = readRDS("DH_seurat.rds")

pdf("DH.umap_plot.pdf")
DimPlot(object = pbmc, reduction = "umap",label=T)
dev.off()

png("DH.violin_plot.png",height=1024,width=2048)
VlnPlot(object = pbmc, features = c("Mog","Apoe", "Pdgfra","C1qb",
  "Slc17a7","Gad2", "Pvalb","Sst","Npy","Vip","Nov","Spink8","Lpl", "Prox1"))
dev.off()


tab = data.frame(sample=pbmc$sample,cluster=pbmc$seurat_clusters)
summ = plyr::count(tab[,1])
tab2 = plyr::count(tab)
tab2$total = summ$freq[match(tab2$sample,summ$x)]
tab2$frac = tab2$freq/tab2$total

pdf("DH_celltype_frac.pdf",height=5,width=10)
ggplot(tab2) + geom_col(aes(x=cluster,y=frac,fill=sample),position="dodge") +
  theme_bw()
dev.off()

pbmc$seurat_clusters.stage = paste0(pbmc$seurat_clusters,"_",pbmc$stage)
Idents(object = pbmc) <- "seurat_clusters.stage"

f = FindMarkers(object = pbmc, ident.1 = "1_03", ident.2 = "1_18", 
    verbose = FALSE)

head(x = f, n = 15)

Idents(pbmc) = "seurat_clusters"
ca1 <- subset(x = pbmc, idents = 1)
Idents(object = ca1) <- "stage"
avg.ca1 <- log1p(x = AverageExpression(object = ca1, verbose = FALSE)$RNA)

p1 <- ggplot(avg.ca1) + geom_point(aes(x=`03`,y=`18`))
p1 <- LabelPoints(plot = p1, points = rownames(f), repel = TRUE)

pdf("CA1.pdf")
p1
dev.off()

FeaturePlot(object = pbmc, features="AC149090.1")


