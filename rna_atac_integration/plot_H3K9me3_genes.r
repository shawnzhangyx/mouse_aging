library(Seurat)
library(ggplot2)

setwd("../../analysis/rna_atac_integration")

pbmc = readRDS("/projects/ps-renlab/lamaral/projects/aging_RNA/DH/analysis/DH_seurat_rmdoub_filtered.rds")
meta =read.table("../../scripts/rna_atac_integration/rna_cell_type.consistent.txt",header=T)
#pbmc$seurat_clusters.stage = paste0(pbmc$seurat_clusters,"_",pbmc$stage)
#Idents(object = pbmc) <- "seurat_clusters.stage"
pbmc$celltype = meta$celltype[match(pbmc$seurat_clusters,meta$cluster)]

pbmc$celltype.stage = paste0(pbmc$celltype,"_",pbmc$stage)
sequence = sort(unique(pbmc$celltype.stage))
pbmc$celltype.stage = factor(pbmc$celltype.stage, levels=sequence)
Idents(pbmc) = "celltype.stage"
pbmc = pbmc[,!is.na(pbmc$celltype)]


dot = DotPlot(pbmc,features=c("Setdb1","Setdb2","Ehmt1","Ehmt2","Suv39h1","Suv39h2","Kdm3a","Kdm3b","Kdm4a","Kdm4b","Kdm4c","Kdm4d"))

dat = dot$data
dat$id = factor(dat$id, levels=sequence)
dat$ct = sub("(.*)_(..)","\\1",dat$id)
dat$class = "Glia"
dat$class[which(dat$ct %in% c("DG","CA1","CA23","Sub_Ent"))] = "ExN"
dat$class[which(dat$ct %in% c("Inh"))] = "InN"

g1 =  ggplot(dat) + geom_col(aes(x=id,y=avg.exp,fill=ct)) + 
  facet_grid(features.plot~class,scale="free",space="free_x",) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

g2 = ggplot(dat) + geom_col(aes(x=id,y=pct.exp,fill=ct)) +
  facet_grid(features.plot~class,scale="free",space="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("DH.H3K9me3_levels.by_celltype_stage.pdf",height=12,width=10)
g1
g2
dev.off()
#VlnPlot(pbmc,features=c("Setdb1"))
#prdms = rownames(pbmc)[grep("Prdm",rownames(pbmc))]
#dot = DotPlot(pbmc,features=prdms)
#cbx = #rownames(pbmc)[grep("Cbx",rownames(pbmc))]
#dot = DotPlot(pbmc,features=cbx)


pbmc = readRDS("/projects/ps-renlab/lamaral/projects/aging_RNA/FC/analysis/FC_seurat_rmdoub_filtered.rds")
meta = read.csv("/projects/ps-renlab/lamaral/projects/aging_RNA/FC/analysis/rna_atac_cons_key.csv", header =F)

#meta =read.table("../../scripts/rna_atac_integration/rna_cell_type.consistent.txt",header=T)
#pbmc$seurat_clusters.stage = paste0(pbmc$seurat_clusters,"_",pbmc$stage)
#Idents(object = pbmc) <- "seurat_clusters.stage"
pbmc$celltype = meta$V2[match(pbmc$RNA_snn_res.0.3,meta$V1)]

pbmc$celltype.stage = paste0(pbmc$celltype,"_",pbmc$stage)
sequence = sort(unique(pbmc$celltype.stage))
pbmc$celltype.stage = factor(pbmc$celltype.stage, levels=sequence)
Idents(pbmc) = "celltype.stage"
pbmc = pbmc[,!is.na(pbmc$celltype)]
dot = DotPlot(pbmc,features=c("Setdb1","Setdb2","Ehmt1","Ehmt2","Suv39h1","Suv39h2","Kdm3a","Kdm3b","Kdm4a","Kdm4b","Kdm4c","Kdm4d"))
dat = dot$data

dat$id = factor(dat$id, levels=sequence)
dat$ct = sub("(.*)_(..)","\\1",dat$id)
dat$class = "Glia"
dat$class[which(dat$ct %in% c("Claustrum","L2.3","L4","L5.Deptor","L5.Fezf2","L5.Parm1","L6"))] = "ExN"
dat$class[grep("InN", dat$ct)] = "InN"

g1 =  ggplot(dat) + geom_col(aes(x=id,y=avg.exp,fill=ct)) + 
  facet_grid(features.plot~class,scale="free",space="free_x",) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

g2 = ggplot(dat) + geom_col(aes(x=id,y=pct.exp,fill=ct)) +
  facet_grid(features.plot~class,scale="free",space="free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("FC.H3K9me3_levels.by_celltype_stage.pdf",height=12,width=10)
g1
g2
dev.off()

