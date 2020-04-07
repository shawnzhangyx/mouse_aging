tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]
setwd(paste0("../../analysis/snapATAC/",tissue))
system("mkdir cluster_feature")

a=data.frame(fread(paste0(tissue,".peaks.counts"),skip=1))
b =read.delim(paste0(tissue,".peaks.counts.summary"))
cnts = a[,-c(1:6)]

grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..sorted.bam","\\1",colnames(cnts)))

library(edgeR)
library(doParallel)
registerDoParallel(cores=rank)

foreach(grp = unique(grps)) %dopar% {
print(grp)
test = rowSums(cnts[,which(grps == grp)])
ctrl = rowSums(cnts[,which(grps != grp)])
dat = cbind(ctrl,test)
rownames(dat) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)

group <- factor(c(1,2));
design <- model.matrix(~group);
y <- DGEList(counts=dat, group=group);
y = calcNormFactors(y)
bcv = 0.1
ext = exactTest(y, dispersion=bcv^2)

fdr = p.adjust(ext$table$PValue,method="BH")
out = cbind(ext$table,fdr)
out = out[order(out$logFC<0,out$PValue),]

write.table(out,paste0("cluster_feature/",tissue,".",grp,".peaks.sortby_logFC.txt"),sep="\t",quote=F)
}

## only look at the promoters
a=read.delim(paste0(tissue,".promoters.counts"),skip=1)
b =read.delim(paste0(tissue,".promoters.counts.summary"))
cnts = a[,-c(1:6)]

grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..sorted.bam","\\1",colnames(cnts)))

foreach(grp=unique(grps)) %dopar% {
print(grp)

test = rowSums(cnts[,which(grps == grp)])
ctrl = rowSums(cnts[,which(grps != grp)])
dat = cbind(ctrl,test)
rownames(dat) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)

group <- factor(c(1,2));
design <- model.matrix(~group);
y <- DGEList(counts=dat, group=group);
y = calcNormFactors(y)
bcv = 0.1
ext = exactTest(y, dispersion=bcv^2)

fdr = p.adjust(ext$table$PValue,method="BH")
out = cbind(ext$table,fdr)
out = out[order(out$logFC<0,out$PValue),]

write.table(out,paste0("cluster_feature/",tissue,".",grp,".prom.sortby_logFC.txt"),sep="\t",quote=F)
}

