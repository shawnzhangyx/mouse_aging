a=read.delim("luisa_all_tissues/all_merged_40kL.cluster.meta.txt")
m1 = read.csv("cluster_keys.csv")
m2 = read.delim("../celltype_annotation.txt")
m2 = m2[!is.na(m2$Clade),]


a$ct_all = m1$ct[match(a$cluster,m1$cluster)]
a$ct_tissue = paste0(m2$Tissue,".",m2$Name)[match(a$tissue_cluster,paste0(m2$Tissue,"_",m2$Cluster))]


library(SnapATAC)

load("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/ATAC_imputed.RData")

dh = data.frame(sample=p.x.sp@sample, barcode=p.x.sp@barcode, rna=p.x.sp@cluster,score=p.x.sp@metaData$predict.max.score,TSSenrich=p.x.sp@metaData$TSS_enrich,TotalReads=p.x.sp@metaData$UQ) 

dh$ct_all = a$ct_all[match(paste(dh$sample,dh$barcode),paste(a$sample,a$barcode))]
dh$ct_tissue = a$ct_tissue[match(paste(dh$sample,dh$barcode),paste(a$sample,a$barcode))]
m3 = read.table("../../../scripts/rna_atac_integration/rna_cell_type.txt",header=T)
dh$ct_rna = m3$celltype[match(dh$rna,m3$cluster)]

table(dh$ct_rna,dh$ct_tissue)
table(dh$ct_rna,dh$ct_all)
table(dh$ct_tissue,dh$ct_all)

dh$all_tissue_match="no"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Asc" & dh$ct_all=="Astrocyte")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.CA1" & dh$ct_all=="CA")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.CA2.3" & dh$ct_all=="CA")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.DG.1" & dh$ct_all=="DG")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.DG.2" & dh$ct_all=="DG")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.DG.3" & dh$ct_all=="DG")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Endo" & dh$ct_all=="Endothelial")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Inh" & dh$ct_all=="InN")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Mgc" & dh$ct_all=="micglia")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Ogc" & dh$ct_all=="Oligodendrocytes")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Opc" & dh$ct_all=="Opc")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Peri" & dh$ct_all=="Pericyte")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.SMC" & dh$ct_all=="SmoothMuscle")]="yes"
dh$all_tissue_match[which(dh$ct_tissue=="DH.Sub_Ent" & dh$ct_all %in% c("L2.3.4","L5","L6"))]="yes"


#mis = dh[which(dh$ct_tissue=="DH.Asc" & dh$ct_all=="DG"),]
#ggplot(dh) + geom_violin(aes(x=all_tissue_match,y=score))

#ct rna  tissue all
#        Asc  1861  1811  1786
#        CA1  4266  4201  4233
#     CA23.1  1428  1428  1406
#     CA23.2   341   324   319
#         DG 14847 13585 14648
#       Endo   486   483   480
#        Inh  1867  1663  1645
#        Mgc  1447  1444  1403
#        Ogc  7226  7215  7146
#        Opc   887   881   872
#       Peri   436   435   381
#        SMC   214   213   205
#  Sub_Ent.1   899   899   813
#  Sub_Ent.2   499   494   447
#  Sub_Ent.3   486   482   439
## write the table into tmp.rna.tissue.all.txt
tmp = read.delim("tmp.rna.tissue.all.txt")
tmp$tissue_rate = tmp$tissue/tmp$rna
tmp$all_rate = tmp$all/tmp$rna

melted = melt(tmp[,c("ct","tissue_rate","all_rate")])


g1 = ggplot(melted) + geom_col(aes(x=ct,y=value,fill=variable),position="dodge")+
  theme_bw(base_size=15)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


g2 = ggplot(dh) + geom_boxplot(aes(x=all_tissue_match,y=score),outlier.shape=NA) +
  theme_bw(base_size=20)
g3 = ggplot(dh) + geom_boxplot(aes(x=all_tissue_match,y=TSSenrich),outlier.shape=NA) +
  theme_bw(base_size=20)
g4 = ggplot(dh) + geom_boxplot(aes(x=all_tissue_match,y=TotalReads),outlier.shape=NA) +
  theme_bw(base_size=20)


pdf("DH_celltype_classification_by_tissue_all.vs.rna_label_transfer.pdf")
g1
g2
g3
g4
dev.off()


