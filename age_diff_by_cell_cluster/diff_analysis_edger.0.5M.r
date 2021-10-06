setwd("../../analysis/snapATAC/de_peaks.subsample")
library(edgeR)

reads = read.delim("select_celltypes.0.5M.subsample.counts",skip=1)

cnts = reads[,-c(1:6)]

ct = sub("bam.subsample.0.5M.(.*).(..).rep..bam","\\1",colnames(cnts))
age = sub("bam.subsample.0.5M.(.*).(..).rep..bam","\\2",colnames(cnts))
## unique cell types)
cts = unique(ct)

system("mkdir age_diff_edgeR.0.5M")
for (ict in cts) {
print(ict)
counts = cnts[,which(ct == ict)]

grps = age[which(ct ==ict)]

y = DGEList(counts)
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)
design = model.matrix(~0+grps)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)


lrt = glmLRT(fit_tag, contrast =c(-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]

write.csv(out, paste0("age_diff_edgeR.0.5M/",ict,".edger.txt"))
}

sig_list = NULL
for (ict in cts) {
tmp = read.csv(paste0("age_diff_edgeR.0.5M/",ict,".edger.txt"))
sig = c( length(which(tmp$PValue<0.01)), length(which(tmp$PValue<0.001)),length(which(tmp$PValue<0.0001)))
sig_list = c(sig_list,sig)
}

dat = data.frame(ct=rep(cts,each=3),de=sig_list,pval=rep(c(0.01,0.001,0.0001),rep=length(cts)),stringsAsFactors=F)

meta =read.delim("../../../aging_share/figures/celltype_annotation.txt")
dat$name = paste0(meta$Tissue,".",meta$Name)[match(dat$ct, paste0(meta$Tissue,".",meta$Cluster))]
dat$class = meta$Clade[match(dat$ct, paste0(meta$Tissue,".",meta$Cluster))]
dat$tissue =meta$Tissue[match(dat$ct, paste0(meta$Tissue,".",meta$Cluster))] 
dat$mitotic = meta$Mitotic[match(dat$ct, paste0(meta$Tissue,".",meta$Cluster))]
dat = subset(dat, !is.na(class))
dat = dat[order(dat$de),]
#dat$name = factor(dat$name,levels=unique(dat$name))

library(gridExtra)
g1 = ggplot(dat) + geom_boxplot(aes(x=class,y=de)) + 
  geom_jitter(aes(x=class,y=de),width=0.2) + 
  facet_grid(pval~.,scales="free")  
g2 = ggplot(dat) + geom_boxplot(aes(x=tissue,y=de)) + 
  geom_jitter(aes(x=tissue,y=de),width=0.2) + 
  facet_grid(pval~.,scales="free")  
g3 = ggplot(dat) + geom_boxplot(aes(x=mitotic,y=de)) +
  geom_jitter(aes(x=mitotic,y=de),width=0.2) +
  facet_grid(pval~.,scales="free")

pdf("de_peaks_by_class_tissue.0.5M.pdf",height=8,width=12)
#grid.arrange(grobs=list(g1,g2,g3),ncol=3,layout_matrix=rbind(c(rep(1,8),rep(2,6),rep(3,3))))
grid.arrange(grobs=list(g1,g2,g3),ncol=3,widths=c(8,6,3))
dev.off()

dat1 = subset(dat,pval==0.01)
dat1 = dat1[order(dat1$de),]
dat1$name = factor(dat1$name,levels=dat1$name)
dat2 = subset(dat,pval==0.001)
dat2 = dat2[order(dat2$de),]
dat2$name = factor(dat2$name,levels=dat2$name)
dat3 = subset(dat,pval==0.0001)
dat3 = dat3[order(dat3$de),]
dat3$name = factor(dat3$name,levels=dat3$name)

g4 = ggplot(dat1) + geom_point(aes(x=name,y=de)) + ylim(0,NA) + coord_flip()
g5 = ggplot(dat2) + geom_point(aes(x=name,y=de)) + ylim(0,NA) + coord_flip() 
g6 = ggplot(dat3) + geom_point(aes(x=name,y=de)) + ylim(0,NA) + coord_flip() 
pdf("de_peaks_by_celltypes.0.5M.pdf",height=6,width=12)
grid.arrange(grobs=list(g4,g5,g6), ncol=3)
dev.off()

#ggplot(dat) + geom_col(aes(x=name,y=de,fill=class)) + coord_flip()





