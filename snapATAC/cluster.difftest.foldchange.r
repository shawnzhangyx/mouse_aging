library(SnapATAC)

setwd("../../analysis/snapATAC/DH/snapFiles/")

system("mkdir ../subsample.age_diff_foldchange")
load("DH.pool.snapATAC.cluster.RData")

x.sp = addPmatToSnap(x.sp)

pmat = x.sp@pmat


peak = x.sp@peak$name

cluster = x.sp@cluster
sample = x.sp@sample
cluster_sample =paste(cluster,sample,sep=".")

#out = colSums(pmat[which(cluster_sample=="1.DH_03_rep1"),])
# number of cells to subsample from. 
NUM =200
rank = 15

#cl = "6"

library(doParallel)
registerDoParallel(cores=rank)

foreach(cl=1:15) %dopar% {
#for (cl in 1:15){
print(cl)
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
fc = log( (a18+1)/(a03+1) )
out = cbind(a03,a10,a18,rowSum,fc)

write.csv(out, paste0("../subsample.age_diff_foldchange/",cl,".foldchange.",NUM,".txt"))

}

################### count the number of differential peaks using the same cut-off. 

dat = list()

dat = foreach(cl=1:rank) %dopar% {
tmp = read.csv(paste0("../subsample.age_diff_foldchange/",cl,".foldchange.",NUM,".txt"))
tmp$cluster = cl
tmp
}


dat2 = do.call(rbind, dat)

ggplot(subset(dat2,rowSum<100)) + geom_histogram(aes(x=rowSum))

t2 = table(dat2$rowSum>=20 & abs(dat2$fc)>=log(2), dat2$cluster)
mt2 =melt(t2)

ncell = table(cluster)
mt2$ncells = ncell[match(mt2$Var2,names(ncell))]

pdf("../subsample.age_diff_foldchange/number_of_diff_peaks.pdf")
ggplot(subset(mt2,Var1==TRUE)) + geom_col(aes(x=Var2,y=value))
ggplot(subset(mt2,Var1==TRUE)) + geom_text(aes(y=value,x=ncells,label=Var2)) + 
  xlim(0,10000) + ylim(0,16000) + xlab("Number of cells") + 
  ylab("Number of Differential Peaks") 

dev.off()




