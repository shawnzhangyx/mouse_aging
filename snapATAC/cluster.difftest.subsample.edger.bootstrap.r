library(SnapATAC)
tissue="HT"
rank = 15

setwd(paste0("../../analysis/snapATAC/",tissue,"/snapFiles/"))
load(paste0(tissue,".pool.snapATAC.cluster.RData"))
set.seed(1997)

pmat = x.sp@pmat


peak = x.sp@peak$name

cluster = x.sp@cluster
sample = x.sp@sample
cluster_sample =paste(cluster,sample,sep=".")

#out = colSums(pmat[which(cluster_sample=="1.DH_03_rep1"),])
# number of cells to subsample from. 
NUM =200

library(edgeR)
library(doParallel)
registerDoParallel(cores=10)

num_diff_out = foreach(cl=1:rank) %dopar% {
# repeat 20 times. 
num_diff = NULL
for (i in 1:20){

print(paste(cl,i))
samples = unique(sample)
dat = list()
for (ss in samples) { 
  idx1 = which(cluster_sample==paste(cl,ss,sep="."))
  idx2 = sample(idx1,NUM,replace=T)
  dat[[ss]] = colSums(pmat[idx2,])
  }

mat = do.call(cbind,dat)
rownames(mat) = peak

grps = substr(colnames(mat),4,5)

y = DGEList(mat)
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)
design = model.matrix(~0+grps)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

contrast.matrix = matrix(c(c(-1,1,0),c(0,-1,1),c(-1,0,1)),nrow=3)
lrt = glmLRT(fit_tag, contrast =contrast.matrix)
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]

num_diff = c(num_diff, length(which(out$PValue<0.01)) )
}
num_diff
}


tab = melt(num_diff_out)
ncell = table(cluster)
tab$ncells = ncell[match(tab$L1,names(ncell))]
## total number of reads: 
rsums = rowSums(pmat)
data = data.frame(cluster,rsums)
readmean = aggregate(rsums~cluster,data,mean)
tab$readmean = readmean$rsums[match(tab$L1,readmean$cluster)]


save(tab, "../subsample.age_diff_edger/tab.Rdata")

system("mkdir ../subsample.age_diff_edger")
pdf("../subsample.age_diff_edger/number_of_diff_peaks.bootstrap.pdf")

ggplot(tab) + geom_boxplot(aes(x=factor(L1),y=value))
ggplot(tab) + geom_point(aes(ncells,value,color=factor(L1))) + 
  geom_text(data=subset(tab,!duplicated(L1)),aes(ncells,value,label=L1),size=5) + 
  ylab("Number of Differential Peaks")
ggplot(tab) + geom_point(aes(readmean,value,color=factor(L1))) +
  geom_text(data=subset(tab,!duplicated(L1)),aes(readmean,value,label=L1),size=5) +
  ylab("Number of Differential Peaks")
  
dev.off()




