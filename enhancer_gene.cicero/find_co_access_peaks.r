tissue=commandArgs(trailing=T)[1]
setwd("../../analysis/cicero_results/")
pr = read.delim(paste0(tissue,"/",tissue,"_peak_overlap_promoter.txt"),header=F,stringsAsFactors=F)
co = read.csv(paste0(tissue,"/",tissue,".cicero_conns.25.csv"),stringsAsFactors=F)

pr$name = paste0(pr$V1,"_",pr$V2,"_",pr$V3)
# find co-access peaks.
co2 = co[which(co$Peak1 %in% pr$name),]
co2$gene = pr$V9[match(co2$Peak1,pr$name)]
co2$chr = sub("(.*)_(.*)_(.*)","\\1",co2$Peak2)
co2$start = sub("(.*)_(.*)_(.*)","\\2",co2$Peak2)
co2$end = sub("(.*)_(.*)_(.*)","\\3",co2$Peak2)

write.table(co2[,c("chr","start","end","Peak1","gene")],paste0(tissue,"/",tissue,"_peak_coaccess_toGene.25.bed"),row.names=F,col.names=F,sep="\t",quote=F)


