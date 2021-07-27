setwd("../../analysis/rna_atac_integration")

a=read.delim("concordant_changes_in_rna_atac.500kb.txt")

b = read.delim("concordant_changes_in_rna_atac.txt")

o = a[-which(paste(a$peak,a$gene,a$ct) %in% paste(b$peak,b$gene,b$ct)),]


a = a[order(a$cor.fdr),]
a = a[order(-a$atac.logP),]
a = a[order(-a$rna.p.adj),]

# remove duplicated gene-cell type pairs. 
b = a[!duplicated(a[,c(2,4)]),]

## tail(sort(table(b$gene)),30)

a[which(a$gene=="Ptgds"),]
# Nrxn3
## Cd9 in Mgc.

#table(table(a$peak))
#
#  1   2   3   4   6   7   8
#  439  22   2   3   1   1   1

# a[which(a$peak=="chr13 3357756 3358757"),]


d = a[!duplicated(a$gene),]


genage = read.delim("../../annotations/genage.mouse.symbol.txt",stringsAsFactors=F)

g = a[which(a$gene %in% genage$mt.Co1),]

# Cd9 Mgc increase  #peak not good

# Robo1 DG decrease 
# Apoe Ogc increase
# Nrg1 CA1 decrease DG increase

a$genage = a$gene %in% genage$mt.Co1

a[which(a$gene == "Robo1"),]
#a[which(a$gene == "Apoe"),]
a[which(a$gene == "Nrg1"),]


