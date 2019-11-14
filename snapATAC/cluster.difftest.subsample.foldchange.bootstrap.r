library(SnapATAC)

setwd("../../analysis/snapATAC/DH/snapFiles/")
load("DH.pool.snapATAC.cluster.RData")
set.seed(1997)

pmat = x.sp@pmat


peak = x.sp@peak$name

cluster = x.sp@cluster
sample = x.sp@sample
cluster_sample =paste(cluster,sample,sep=".")

#out = colSums(pmat[which(cluster_sample=="1.DH_03_rep1"),])
# number of cells to subsample from. 
NUM =200

library(doParallel)
registerDoParallel(cores=10)

num_diff_out = foreach(cl=1:15) %dopar% {
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

a03 = rowSums(mat[,grps=="03"])
a10 = rowSums(mat[,grps=="10"])
a18 = rowSums(mat[,grps=="18"])
rowSum = rowSums(mat)
fc = log( (a18/sum(a18)*1e6+1)/(a03/sum(a03)*1e6+1) )
out = data.frame(a03,a10,a18,rowSum,fc)


num_diff = c(num_diff, length(which(out$rowSum>20 & abs(out$fc)>=log(2) )) )
}
num_diff
}


tab = melt(num_diff_out)
ncell = table(cluster)
tab$ncells = ncell[match(tab$L1,names(ncell))]

pdf("../subsample.age_diff_foldchange/number_of_diff_peaks.bootstrap.pdf")

ggplot(tab) + geom_boxplot(aes(x=factor(L1),y=value))
ggplot(tab) + geom_point(aes(ncells,value,color=factor(L1))) + 
  geom_text(data=subset(tab,!duplicated(L1)),aes(ncells,value,label=L1),size=5) + 
  ylab("Number of Differential Peaks")
dev.off()




