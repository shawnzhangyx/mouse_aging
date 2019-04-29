setwd("../../analysis/Yang_NMF_method/heart")

a = read.delim("../../../analysis/Yang_NMF_method/heart/R10/heart.R10.statH")
b= read.table("../../../data/snATAC/ct2peaks/heart.all.barcode.peak.cnts")
a$ripeak = b$V2[match(a$X..xgi,b$V1)]

b2 = read.table("../../../data/snATAC/bam_bowtie2_Olivier/qc/heart.all.dedup.filter.barcode.cnts")
a$rdf = b2$V2[match(a$X..xgi,b2$V1)]
  
#ggplot(a) + geom_violin(aes(x=factor(class0+1),y=ripeak)) + scale_y_log10()

pdf("diagnostic.cluster.reads.pdf")
ggplot(a) + geom_boxplot(aes(x=factor(class0+1),y=ripeak)) + scale_y_log10()

ggplot(a) + geom_boxplot(aes(x=factor(class0+1),y=rdf)) + scale_y_log10()

ggplot(a) + geom_boxplot(aes(x=factor(class0+1),y=ripeak/rdf)) + scale_y_log10()
dev.off()