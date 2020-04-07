library(ggrepel)
library(gridExtra)
library(edgeR)

tissue=commandArgs(trailing=T)[1]
tissue="FC"

setwd("../../analysis/repeats_RepEnrich2/RE_out/")
## read meta info
meta = list.files(path="../metaInfo",recursive=T,full.names=T,pattern="txt")
dat_list = list()
for (f in meta) {
dat_list[[f]] = read.table(f)
}

meta = do.call(rbind,dat_list)
meta$sample = paste0(meta$V1,".",meta$V9)
total = aggregate(V3~sample,meta,sum)
total = total[grep(tissue,total$sample),]

## read the matrix
fs = list.files(path=".",recursive=T,pattern="[0-9]_fraction_counts.txt")
fs = fs[grep(tissue,fs)]
dat_list = list()
for (f in fs) {
dat_list[[f]] = read.table(f)
}



dat = Reduce(function(...) merge(...,by = c("V1","V2","V3")),
       dat_list)
sample = sub("(.*)\\/(.*)","\\1",fs)
grps = sub(".._(..)_rep..(.*)","\\1",sample)
cluster = sub(".._(..)_rep..(.*)","\\2",sample)

colnames(dat) = c("family","class","name",sample)
rownames(dat) = dat$family

for (clu in sort(unique(cluster)) ) {

idx = which(cluster == clu)+3
y = DGEList(counts=dat[,idx])
y$samples$lib.size = total[idx-3,"V3"]
#y = calcNormFactors(y)
groups = grps[idx-3] # c("03","03","10","10","18","18")
design = model.matrix(~0+groups)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast =c(-1,0,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]
#out = out[-which(rownames(out) %in% c("Others")),]
print(c(clu, length(which(out$fdr<0.01))))
}


ggplot(out) + geom_point(aes(x=logCPM,y=logFC,color=fdr<0.01))

ggplot(out) + geom_point(aes(x=logFC,y=-log10(fdr),color=fdr<0.01))



