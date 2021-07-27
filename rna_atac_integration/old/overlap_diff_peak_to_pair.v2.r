link=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/weighted_regression.txt"))

setwd("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/snapATAC/DH/age_diff_edgeR.snap/")

files = list.files(pattern="*(up|down).bed")

dat.list = list()

for (file in files) { 
  dat.list[[file]] = data.frame(fread(file))
  dat.list[[file]]$file = file
  }

atac = do.call(rbind,dat.list)
atac$peak = paste(atac$V1,atac$V2,atac$V3)
atac$cluster = sub("(.*)\\.(.*)\\.bed","\\1",atac$file)
atac$direction = sub("(.*)\\.(.*)\\.bed","\\2",atac$file)

atacm = read.delim("../../../../aging_share/figures/celltype_annotation.txt")
atacm = atacm[which(atacm$Tissue=="DH"),]
atacm$Name = c("Ogc","DG","DG","CA1","DG","Inh","Sub_Ent","Asc","CA23","Mgc","CP","Opc","Endo","Peri","SMC")
atacm = atacm[-11,]
atac$ct = atacm$Name[match(atac$cluster,atacm$Cluster)]

link$peak = paste(link$ATAC_Chr,link$ATAC_start,link$ATAC_end)
link$fdr = p.adjust(link$Pvalue,method="BH")
# only limit to peaks withfdr < 0.05
link = link[link$fdr<0.05,]

## overlap with diff peaks. 
b = link[which(link$peak %in% atac$peak),]
peaks = atac[which(atac$peak %in% link$peak),]
peaks$gene = b$Geneid[match(peaks$peak,b$peak)] # there is possibility that one peak get assigned to multiple gene. 

setwd("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/rna_atac_integration")


fs2 = list.files(pattern="3vs18_nocutoff.txt",path="RNA_diff",full.names=T)

dat.list = list()

for (file in fs2) {
  dat.list[[file]] = data.frame(fread(file))
  dat.list[[file]]$ct = sub(".*\\/(.*)_03vs18_nocutoff.txt","\\1",file)
    }
rna = do.call(rbind,dat.list)
rownames(rna) = NULL
rna = rna[which(rna$p_val_adj<0.05),]

rna$p_val = -log10(rna$p_val)
rna$p_val_adj = -log10(rna$p_val_adj)
rna$direction = ifelse(rna$avg_logFC<0,"up","down")
colnames(rna)[1] = "Geneid"
genes = unique(rna$Geneid)

pbmc = readRDS("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/DH_seurat_rmdoub_filtered.rds")
meta =read.table("../../scripts/rna_atac_integration/rna_cell_type.consistent.txt",header=T)
pbmc$celltype = meta$celltype[match(pbmc$seurat_clusters,meta$cluster)]
pbmc =pbmc[,!is.na(pbmc$celltype)]
Idents(object = pbmc) <- "celltype"

dot = DotPlot(pbmc,features=genes)

dot.data = dot$data

rna$avgRNA = dot.data$avg.exp[match(paste(rna$Geneid,rna$ct), paste(dot.data$features.plot,dot.data$id))]


#> table(rna$ct)
#
#    Asc     CA1    CA23      DG    Endo     Inh     Mgc     Ogc     Opc    Peri
#     14     254      32     414       2     396      17     181       9       1
#    SMC Sub_Ent
#       1     130
#> table(atac$ct)
#    Asc     CA1    CA23      DG    Endo     Inh     Mgc     Ogc     Opc     Peri
#    619    2551    2101    7154      81    1693     484    1746     469     200
#    SMC Sub_Ent
#     45    2209

cts = unique(rna$ct)
## for each cell type find the target genes of DE peaks. and then ask how many genes are differential. 

ct = "DG"
dir = "down"

out_list = list()
concord_list = list()
 for (ct in cts) {
  for (dir in c("down","up")) {

  peaks = atac[which(atac$ct== ct & atac$direction==dir),]
  tmp.links = link[which(link$peak %in% peaks$peak),]
  genes = unique(tmp.links$Geneid)  
  tmp.rna = rna[which(rna$Geneid %in% genes),]
  ct.rna = tmp.rna[which(tmp.rna$ct ==ct),]
  down = length(which(ct.rna$direction=="down"))
  up = length(which(ct.rna$direction=="up"))
  combined =  c(ct,dir,length(genes),down,up)
  out_list[[length(out_list)+1]] = combined

  concord.rna = tmp.rna[which(tmp.rna$ct ==ct & tmp.rna$direction==dir),]
  concord.links = tmp.links[which(tmp.links$Geneid %in% concord.rna$Geneid),]
  if(nrow(concord.links) ==0) { next() }
  concord.links$GenePos = paste0(concord.links$Chr,":",concord.links$Start,"-",concord.links$End)
  concord.rna.links = merge(concord.rna[,-2],concord.links[,-c(2:11)],by="Geneid")
  concord.atac = peaks[which(peaks$peak %in% concord.rna.links$peak),c(4,5,7)]
  concord.rna.links.atac = merge(concord.rna.links,concord.atac,by="peak")
  concord_list[[length(concord_list)+1]] = concord.rna.links.atac
}
}


tab = data.frame(do.call(rbind,out_list),stringsAsFactors=F)
colnames(tab) = c("ct","dirATAC","Ngenes","Ndown","Nup")
tab$Ngenes = as.numeric(tab$Ngenes)
tab$Ndown = as.numeric(tab$Ndown)
tab$Nup = as.numeric(tab$Nup)


Nd = sum(tab$Ngenes[which(tab$dirATAC=="down")])
Nu = sum(tab$Ngenes[which(tab$dirATAC=="up")])


dd = sum(tab$Ndown[which(tab$dirATAC=="down")])
du = sum(tab$Nup[which(tab$dirATAC=="down")])
ud = sum(tab$Ndown[which(tab$dirATAC=="up")])
uu = sum(tab$Nup[which(tab$dirATAC=="up")])

mat = matrix(c(dd,du,ud,uu),nrow=2)
mat
# col: change in ATAC. row: change in RNA. (down,up)
#     [,1] [,2]
#     [1,]  114  105
#     [2,]   81  190
fisher.test(mat)
# p-value = 0.0000007856

concord.out = do.call(rbind,concord_list)
colnames(concord.out) = c("peak","gene","rna.fc","rna.pct1","rna.pct2","rna.p.adj","ct","direction","rna.avg","cor.fdr","GenePos","atac.fc","atac.logP")
concord.out = concord.out[,c(1,2,11,7,8,9,3:6,10,12,13)]
concord.out$rna.fc = -concord.out$rna.fc

write.table(concord.out,"concordant_changes_in_rna_atac.txt",row.names=F,quote=F,sep="\t")



