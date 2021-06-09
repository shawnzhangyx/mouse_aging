library(SnapATAC)

tissue=commandArgs(trailing=T)[1] #"DH"
RData =commandArgs(trailing=T)[2] #"DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData"

setwd(paste0("../../analysis/snapATAC/",tissue,"/snapFiles/"))
#load("DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.cluster.RData")
load(RData)
system("mkdir -p ../age_diff_edgeR.snap/")

pmat = x.after.sp@pmat
peak = x.after.sp@peak$name
cluster = x.after.sp@cluster
sample = x.after.sp@sample
cluster_sample =paste(cluster,sample,sep=".")

#cl = "6"
max_cluster = max(as.numeric(cluster))
library(edgeR)
library(doParallel)
registerDoParallel(cores=max_cluster)

foreach(cl=1:max_cluster) %dopar% {
print(cl)
samples = unique(sample)
dat = list()
for (ss in samples) { 
  idx1 = which(cluster_sample==paste(cl,ss,sep="."))
  dat[[ss]] = colSums(pmat[idx1,])
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

#contrast.matrix = matrix(c(c(-1,1,0),c(0,-1,1),c(-1,0,1)),nrow=3)
#lrt = glmLRT(fit_tag, contrast =contrast.matrix)
lrt = glmLRT(fit_tag, contrast =c(-1,0,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]

write.csv(out, paste0("../age_diff_edgeR.snap/",cl,".edger.txt"))
}

## write differential peaks into bed files. 
for (cl in 1:max_cluster) {
  print(cl)
  tmp = read.csv(paste0("../age_diff_edgeR.snap/",cl,".edger.txt"))
  sig = tmp[which(tmp$PValue<0.01),]
#  sig = tmp[which(tmp$PValue< quantile(tmp$PValue,0.01)),]
  chr= sub("(chr.*):(.*)-(.*)","\\1",sig$X)
  start = sub("(chr.*):(.*)-(.*)","\\2",sig$X)
  end = sub("(chr.*):(.*)-(.*)","\\3",sig$X)
  out = data.frame(chr,start,end,sig$logFC,as.integer(-log10(sig$PValue)))
#  write.table(out, paste0("../age_diff_edgeR.snap/",cl,".both.bed"),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(subset(out,sig.logFC>0), paste0("../age_diff_edgeR.snap/",cl,".up.bed"),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(subset(out,sig.logFC<0), paste0("../age_diff_edgeR.snap/",cl,".down.bed"),row.names=F,col.names=F,quote=F,sep="\t")

  }


