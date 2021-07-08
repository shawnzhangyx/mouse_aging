a=read.delim("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/Robo1_corr_table.txt")
#a$Peakid=NULL
a$Chr=a$Start=a$End=NULL
a = a[order(a$Pvalue),]

# peaks/DH_peak_83639    chr16   72924125 72925126    #0.6052693
b=read.table("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/final_pairs_table_concat.txt",header=T)
b = b[which(b$gene=="Robo1"),]
b = b[order(-b$log10Pval),]
# chr16 72810662 72811663 #170.87230


rna = read.delim("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/rna_counts_consistent.txt")

rna[,-c(1:6)] = sweep(rna[,-c(1:6)],2,colSums(rna[,-c(1:6)]),"/")
rrna = melt(rna[which(rna$Geneid=="Robo1"),-c(1:6)])

ggplot(rrna) + geom_text(aes(x=variable,y=value,label=variable))



dna = read.delim("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/atac_counts_consistent.txt")
dna[,-c(1:6)] = sweep(dna[,-c(1:6)],2,colSums(dna[,-c(1:6)]),"/")


library(ggrepel)

ddna = melt(dna[which(dna$Geneid=="peaks/DH_peak_83639"),-c(1:6)])
rddna = merge(rrna,ddna,by="variable")
ggplot(rddna) + geom_point(aes(x=value.x,y=value.y)) +
  geom_text_repel(aes(x=value.x,y=value.y,label=variable))


ddna = melt(dna[which(dna$Start==72810662),-c(1:6)])
rddna = merge(rrna,ddna,by="variable")
ggplot(rddna) + geom_point(aes(x=value.x,y=value.y))+
  geom_text_repel(aes(x=value.x,y=value.y,label=variable))


# chr16 72772036 72773037
ddna = melt(dna[which(dna$Start==72772036),-c(1:6)])
rddna = merge(rrna,ddna,by="variable")
ggplot(rddna) + geom_point(aes(x=value.x,y=value.y))+
  geom_text_repel(aes(x=value.x,y=value.y,label=variable))


ddna = melt(dna[which(dna$Geneid=="peaks/DH_peak_83584"),-c(1:6)])
rddna = merge(rrna,ddna,by="variable")
ggplot(rddna) + geom_point(aes(x=value.x,y=value.y)) +
  geom_text_repel(aes(x=value.x,y=value.y,label=variable))

dh = read.delim("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/snapATAC/DH/DH.pool.barcode.meta_info.txt")
#meta = read.delim("/projects/ps-renlab/yanxiao/projects/mouse_aging/aging_share/figures/celltype_annotation.txt")
cl_info = data.frame(ct=c("Asc","CA1","CA2.3","DG","DG","DG","Endo","Inh","Mgc","Ogc","Opc","Peri","SMC","Sub_Ent"),cluster=c(8,4,9,2,3,5,13,6,10,1,12,14,15,7))

dh$ct = cl_info$ct[match(dh$cluster,cl_info$cluster)]

dh = dh[!is.na(dh$ct),]
dh$stage = sub("DH_(..)_rep.","\\1",dh$sample)
dh$ct_age = paste0(dh$ct,"_",dh$stage)

agg = aggregate(ct~ct_age,dh,length)
rddna$weight = agg$ct[match(rddna$variable,agg$ct_age)]

summary(lm(value.x~value.y,rddna,weights=rddna$weight))                                 











