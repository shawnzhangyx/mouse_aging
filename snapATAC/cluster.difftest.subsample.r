library(SnapATAC)

setwd("../../analysis/snapATAC/DH/snapFiles/")
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


#cl = "6"

library(edgeR)
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

write.csv(out, paste0("../subsample.age_diff_edger/",cl,".edger.",NUM,".txt"))

}

################### count the number of differential peaks using the same cut-off. 

dat = list()

for (cl in 1:15) {
dat[[cl]] = read.csv(paste0("../subsample.age_diff_edger/",cl,".edger.",NUM,".txt"))
dat[[cl]]$cluster = cl
}


dat2 = do.call(rbind, dat)

t2 = table(dat2$PValue < 0.01, dat2$cluster)
mt2 =melt(t2)

ncell = table(cluster)
mt2$ncells = ncell[match(mt2$Var2,names(ncell))]

pdf("../subsample.age_diff_edger/number_of_diff_peaks.pdf")
ggplot(subset(mt2,Var1==TRUE)) + geom_col(aes(x=Var2,y=value))
ggplot(subset(mt2,Var1==TRUE)) + geom_text(aes(y=value,x=ncells,label=Var2)) + 
  xlim(0,10000) + ylim(0,2500) + xlab("Number of cells") + 
  ylab("Number of Differential Peaks") 

dev.off()




