setwd("../../analysis/snapATAC/de_peaks.subsample")
library(edgeR)

reads = read.delim("select_celltypes.grid.subsample.counts",skip=1)

cnts = reads[,-c(1:6)]

ct = sub("bam.subsample.grids.(.*).(..).rep..(.*).bam","\\1",colnames(cnts))
age = sub("bam.subsample.grids.(.*).(..).rep..(.*).bam","\\2",colnames(cnts))
lib = sub("bam.subsample.grids.(.*).(..).rep..(.*).bam","\\3",colnames(cnts))

## unique cell types)
#cts = unique(ct)
ct.lib = paste0(ct,".",lib)

ct.libs = unique(ct.lib)

system("mkdir age_diff_edgeR.grid ")
for (ict in ct.libs) {
print(ict)
counts = cnts[,which(ct.lib == ict)]

grps = age[which(ct.lib ==ict)]

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

write.csv(out, paste0("age_diff_edgeR.grid/",ict,".edger.txt"))
}

sig_list = NULL
for (ict in ct.libs) {
tmp = read.csv(paste0("age_diff_edgeR.grid/",ict,".edger.txt"))
sig = c(length(which(tmp$fdr<0.05)),length(which(tmp$PValue<0.01)), length(which(tmp$PValue<0.001)),length(which(tmp$PValue<0.0001)))
sig_list = c(sig_list,sig)
}

pval = rep(c(0.05, 0.01,0.001, 0.0001),length(ct.libs))
ct = rep(sub("(....).(.*)","\\1",ct.libs),each=4)
lib = rep(sub("(....).(.*)","\\2",ct.libs),each=4)

dat = data.frame(pval,ct,lib,sig_list)

pdf("diff_peaks.by_subsample.grid.pdf")
ggplot(dat) + geom_line(aes(x=lib,y=sig_list,group=ct,color=ct)) +
  facet_grid(pval~.,scales="free_y")
dev.off()


