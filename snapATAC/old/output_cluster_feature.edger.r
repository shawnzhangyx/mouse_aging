tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]
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
registerDoParallel(cores=rank)

y = DGEList(cnts)
y$sample$lib.size=total
y = calcNormFactors(y)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y)

#for (grp in unique(grps) ){
foreach(grp = unique(grps)) %dopar% {
print(grp)
groups = rep("A",length(grps))
groups[which(grps==grp)]= "B"
design = model.matrix(~0+groups)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast = c(-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(lrt$table,fdr)
out = out[order(out$PValue),]
out = out[order(-out$logFC),]

write.table(out,paste0("cluster_feature/",tissue,".",grp,".peaks.sortby_logFC.txt"),sep="\t",quote=F)
}

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


foreach(grp=unique(grps)) %dopar% {
print(grp)
groups = rep("A",length(grps))
groups[which(grps==grp)]= "B"
design = model.matrix(~0+groups)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast = c(-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(lrt$table,fdr)
out = out[order(out$PValue),]
out = out[order(-out$logFC),]

write.table(out,paste0("cluster_feature/",tissue,".",grp,".prom.sortby_logFC.txt"),sep="\t",quote=F)
}

