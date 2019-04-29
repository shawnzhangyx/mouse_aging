tissue=commandArgs(trailing=T)[1]
setwd("../../data/snATAC/")

a=data.frame(fread(paste0("counts/",tissue,".summits_ex1k.counts")))
require(edgeR)

counts = a[,-c(1:6)]
rownames(counts)= a$Geneid
colnames(counts) = sub("bam.*dedup_bam.(.*).bc.dedup.bam","\\1",colnames(counts))
#colnames(counts) = group
#design = model.matrix(~0+group)

y= DGEList(counts=counts)#,group=group)
keep = which(rowSums(cpm(y)>1)>=2)
y = y[keep,]
y =  calcNormFactors(y)

pdf(paste0("../../analysis/bulk_analysis/",tissue,".MDS_plots.pdf"))
plotMDS(y)
dev.off()


group = sub(".*\\.(..)\\.(rep.)","\\2",colnames(counts))
design = model.matrix(~0+group)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast = c(-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(a[keep,1:4],lrt$table,fdr)
out = out[order(out$fdr),]

#ggplot(out) + geom_point(aes(x=logCPM,y=logFC,color=fdr<0.05))
write.table(out,paste0("../../analysis/bulk_analysis/",tissue,".rep2.vs.rep1.edger.txt"),row.names=F,quote=F,sep="\t")


group = sub(".*\\.(..)\\.(rep.)","\\1",colnames(counts))
design = model.matrix(~0+group)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast = c(-1,1,0))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(a[keep,1:4],lrt$table,fdr)
out = out[order(out$fdr),]
write.table(out,paste0("../../analysis/bulk_analysis/",tissue,".10.vs.03.edger.txt"),row.names=F,quote=F,sep="\t")

lrt = glmLRT(fit_tag, contrast = c(0,-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(a[keep,1:4],lrt$table,fdr)
out = out[order(out$fdr),]
write.table(out,paste0("../../analysis/bulk_analysis/",tissue,".18.vs.10.edger.txt"),row.names=F,quote=F,sep="\t")

lrt = glmLRT(fit_tag, contrast = c(-1,0,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(a[keep,1:4],lrt$table,fdr)
out = out[order(out$fdr),]
write.table(out,paste0("../../analysis/bulk_analysis/",tissue,".18.vs.03.edger.txt"),row.names=F,quote=F,sep="\t")

write.table(out[which(out$fdr<0.05),],paste0("../../analysis/bulk_analysis/",tissue,".18.vs.03.edger.sig.txt"),row.names=F,quote=F,sep="\t")


