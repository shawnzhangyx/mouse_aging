tissue=commandArgs(trailing=T)[1]
setwd("../../data/snATAC/")

a=data.frame(fread(paste0("counts/",tissue,".summits_ex1k.counts")))

a$total = rowSums(a[,-c(1:6)])

a = a[order(-a$total),c(1:6,ncol(a))]
ggplot(a) + geom_histogram(aes(log10(total)),breaks=100)

### remove top 0.5% of the windows. 
a = a[-c(1:(nrow(a)*0.005)),]
a = a[order(as.numeric(sub("peak(.*)","\\1",a$Geneid))),]
write.table(a[,c(2:4,1)],paste0("peaks/",tissue,".all.pooled_summits.ext1k.filter.bed"),row.names=F,sep="\t",quote=F,col.names=F)
