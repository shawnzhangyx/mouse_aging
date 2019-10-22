tissue="HT"  # commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/snapATAC/",tissue))
system("mkdir cluster_feature")
a=data.frame(fread(paste0(tissue,".peaks.counts"),skip=1))
b =read.delim(paste0(tissue,".peaks.counts.summary"))
cnts = a[,-c(1:6)]
total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])

grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..sorted.bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)


library(edgeR)
library(doParallel)

y = DGEList(cnts)
y$sample$lib.size=total
y = calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y)

## only look at the promoters
a=read.delim(paste0(tissue,".promoters.counts"),skip=1)
b =read.delim(paste0(tissue,".promoters.counts.summary"))
cnts = a[,-c(1:6)]
total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])

grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..sorted.bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)

y = DGEList(cnts)
y$sample$lib.size=total
y = calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y)


design = model.matrix(~0+grps)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast = c(0,0,-1,rep(0,6),1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(lrt$table,fdr)
out = out[order(out$PValue),]
#out = out[order(-out$logFC),]
write.table(out,paste0("cluster_feature/C2.vs.C9.prom.txt"),sep="\t",quote=F)



