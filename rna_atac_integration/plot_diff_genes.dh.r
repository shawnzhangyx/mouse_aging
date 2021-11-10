library(Seurat)
library(ggplot2)

setwd("../../analysis/rna_atac_integration")

pbmc = readRDS("/projects/ps-renlab/lamaral/projects/aging_RNA/DH/analysis/DH_seurat_rmdoub_filtered.rds")
meta =read.table("../../scripts/rna_atac_integration/rna_cell_type.consistent.txt",header=T)
#pbmc$seurat_clusters.stage = paste0(pbmc$seurat_clusters,"_",pbmc$stage)
#Idents(object = pbmc) <- "seurat_clusters.stage"
pbmc$celltype = meta$celltype[match(pbmc$seurat_clusters,meta$cluster)]

pbmc$celltype.stage = paste0(pbmc$celltype,"_",pbmc$stage)


a = subset(pbmc, celltype %in% c("DG","CA1","Ogc"))


a$celltype.stage = factor(a$celltype.stage, levels= rev(paste( rep(c("Ogc","DG","CA1"),each=3), rep(c("03","10","18"),3),sep="_")))
a$rep = sub(".*_(rep.)","\\1",a$sample)
a$celltype.stage.rep = factor(paste0(a$celltype.stage,"_",a$rep), levels=rev(paste( rep(c("DG","CA1","Ogc"),each=6), rep(rep(c("03","10","18"),each=2),3),rep(c("rep1","rep2"),9), sep="_")))


Idents(a) = "celltype.stage"

#Idents(a) = "celltype.stage.rep"

dot = 
DotPlot(a,features=c("Robo1","Itgb5","Nrg1"),
 # cols = c("yellow","red")
#  cols=viridis(3)
)

# dot$data$celltype.stage = sub("(.*)_rep.","\\1",dot$data$id)
dot$data$celltype.stage  = dot$data$id 
dot$data$ct = sub("(.*)_..","\\1",dot$data$celltype.stage)


ggplot(dot$data) + geom_col(aes(x=celltype.stage,y=avg.exp,fill=ct)) + 
  facet_grid(features.plot~.,scales="free") #+ 
#  coord_flip() + 


g1 = ggplot(subset(dot$data,features.plot=="Robo1")) + 
  geom_col(aes(x=celltype.stage,y=avg.exp,fill=ct),width=0.8) + 
  coord_flip() + theme_bw()

g2 = ggplot(subset(dot$data,features.plot=="Itgb5")) +
  geom_col(aes(x=celltype.stage,y=avg.exp,fill=ct),width=0.8) +
  coord_flip() + theme_bw()

g3 = ggplot(subset(dot$data,features.plot=="Nrg1")) +
  geom_col(aes(x=celltype.stage,y=avg.exp,fill=ct),width=0.8) +
  coord_flip() + theme_bw()


g4 = ggplot(subset(dot$data,features.plot=="Robo1")) +
  geom_col(aes(x=celltype.stage,y=pct.exp,fill=ct),width=0.8) +
  coord_flip() + theme_bw()

g5 = ggplot(subset(dot$data,features.plot=="Itgb5")) +
  geom_col(aes(x=celltype.stage,y=pct.exp,fill=ct),width=0.8) +
  coord_flip() + theme_bw()

g6 = ggplot(subset(dot$data,features.plot=="Nrg1")) +
  geom_col(aes(x=celltype.stage,y=pct.exp,fill=ct),width=0.8) +
  coord_flip() + theme_bw()


library(gridExtra)
pdf("DiffGenes.plot.pdf",height=8,width=6)
grid.arrange(grobs=list(g1,g4,g2,g5,g3,g6),ncol=2)
dev.off()


