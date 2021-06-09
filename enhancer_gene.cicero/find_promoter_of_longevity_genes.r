setwd("../../analysis/cicero_results/")

a=read.table("../../annotations/genage.mouse.symbol.txt")

b=read.table("../../annotations/gencode.vM10.annotation.gene.tss1k.bed")

d =b[which(b$V4 %in% a$V1),]

write.table(d,"genage.tss1k.bed",row.names=F,col.names=F,quote=F,sep="\t")

