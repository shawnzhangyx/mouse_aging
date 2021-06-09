gene = "Igf1"
tissue="LM"
setwd("../../analysis/cicero_results/")

tss.list = read.table("LM/LM_peak_coaccess_toGene.25.bed",stringsAsFactors=F)
tss = unique(tss.list[tss.list$V5 ==gene,"V4"])
#if (length(tss

conns = read.csv("LM/LM.cicero_conns.25.csv",row.names=1,stringsAsFactors=F)
conns2 = conns[which(conns$Peak1 %in% tss | conns$Peak2 %in% tss),]
library(cicero)
gene_annotation = fread("~/annotations/mm10/gencode.vM10.annotation.gtf")
gene_annotation = data.frame(gene_annotation)
gene_annotation$V9=NULL
gene_annotation = gene_annotation[which(gene_annotation$V3=="transcript"),]

colnames(gene_annotation) = c("chromosome","database","type","start","end","T1","strand","T2")


pdf("Igf1.LM.conns.pdf",width=10,height=5)
plot_connections(conns, "chr10", 87400000,88300000, 
#                gene_model = gene_annotation, 
                coaccess_cutoff = .25, 
                connection_width = .5, 
                collapseTranscripts = "longest" )
dev.off()

peaks = unique(c(conns2$Peak1,conns2$Peak2))
write.table(gsub("_","\t",peaks),"Apoe.DH.peaks.bed",row.names=F,col.names=F,quote=F)


#snp.list = read.table("all_tissue.coA_peak.overlapGWAS.25.filtered.txt")

