setwd("../../analysis/cicero_results")

a=read.delim("all_tissue.coA_peak.overlapGWAS.25.txt",header=F,stringsAsFactors=F)

terms = read.delim("../../data/all_aging_traits_EFO.curated.txt",header=F)
terms$V1 = gsub("[^a-zA-Z]",".",terms$V1)

length(unique(a$V9))
table(a$V5,a$V13)
table(a$V5,a$V9)

a = a[which(a$V9 %in% terms$V1),]
write.table(a,"all_tissue.coA_peak.overlapGWAS.25.filtered.txt",row.names=F,col.names=F,sep="\t",quote=F)

snps = a[!duplicated(a$V10),]
snps = snps[,c("V6","V7","V8","V10")]
write.table(snps,"genAge.coA.SNPs.bed",row.names=F,col.names=F,sep="\t",quote=F)

#table(a$V5,a$V9)
tab = table(a$V5,a$V9)
sort(rowSums(tab>0))
# number of traits 
#  Sp1   Nos3   Mapt   Apoe  Ercc2   Fen1
#  3      4      7      9      9     11


tab2 = table(a$V5,a$V10)
sort(rowSums(tab2>0))
# number of SNPs. 
#  Sp1   Esr1   Hic1  Hspa8   Apoe   Mapt
#  3      4      4      8     12     18

table(a$V5,a$V13)
#       DH FC HT LM
#Apoe   24 24  0  0
#Fen1   22 11  0 22
#Mapt   18 15  0 60
#Ercc2  10  0  0  0
#Nos3    0  0  4  4
#Sp1     5  1  0  0
#Hspa8   9  9  0  0
#Hic1    0  0  1 17
#Esr1    0  0 16  8

## remove duplicated SNPs associated with same trait.
au = a[!duplicated(cbind(a$V5,a$V9,a$V10)),]
tab.trait = table(au$V5,au$V9)
tab.trait[tab.trait>=5] = 5
library(pheatmap)
pdf("Genes.by.Traits.SNPcnts.pdf",height=8,width=8)
pheatmap(tab.trait) # + theme( axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

##
ann = read.table("/projects/ps-renlab/lamaral/projects/Aging/GWAS/bed_LD_snp_files_wtag/hglft_genome_wtag.bed",stringsAsFactors=F)
ann$comb = paste(ann$V5,ann$V6, ann$V4,sep="|")
snps$V4 = ann$comb[match(snps$V4,ann$V5)]
snps = snps[!is.na(snps$V4),]
write.table(snps,"genAge.coA.SNPs.details.bed",row.names=F,col.names=F,sep="\t",quote=F)


