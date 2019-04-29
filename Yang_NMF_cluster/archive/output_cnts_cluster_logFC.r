tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/Yang_NMF_method/",tissue))
a=read.delim(paste0(tissue,".counts"),skip=1)
cnts = a[,-c(1:6)]

grps = paste0("C",sub(".*_(.)....rep..bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = a$Geneid

cnts2 =t( apply(cnts,1,function(x){tapply(x,grps,sum)}) )



library(edgeR)
#design = model.matrix(~0+grps)
y = DGEList(cnts2)
y = calcNormFactors(y)

cpms = cpm(y)
#cpms = y$pseudo.counts
colnames(cpms) = unique(grps)
rownames(cpms) = a$Geneid
cpms = cpms+0.01
logFC = cbind(a[,2:4],log2(sweep(cpms,1,apply(cpms,1,mean),'/')))

for (grp in unique(grps) ){
print(grp)
out = logFC[,c(1:3,which(colnames(logFC)==grp))]
out = out[order(-out[,4]),]

write.table(out,paste0(grp,".logFC.txt"),col.names=F,sep="\t",quote=F)
}


## promoters. 
a=read.delim(paste0(tissue,".prom.counts"),skip=1)
cnts = a[,-c(1:6)]

grps = paste0("C",sub(".*_(.)....rep..bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = a$Geneid

cnts2 =t( apply(cnts,1,function(x){tapply(x,grps,sum)}) )



library(edgeR)
#design = model.matrix(~0+grps)
y = DGEList(cnts2)
y = calcNormFactors(y)

cpms = cpm(y)
#cpms = y$pseudo.counts
colnames(cpms) = unique(grps)
rownames(cpms) = a$Geneid
cpms = cpms+0.01
logFC = cbind(a[,2:4],log2(sweep(cpms,1,apply(cpms,1,mean),'/')))

for (grp in unique(grps) ){
print(grp)
out = logFC[,c(1:3,which(colnames(logFC)==grp))]
out = out[order(-out[,4]),]

write.table(out,paste0(grp,".prom.logFC.txt"),col.names=F,sep="\t",quote=F)
}

