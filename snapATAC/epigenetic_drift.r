library(SnapATAC)

setwd("../../analysis/snapATAC/DH/snapFiles/")
load("DH.pool.snapATAC.cluster.RData")

cluster = 1

bmat = x.sp@bmat

samples = x.sp@sample[x.sp@cluster==cluster & rowSums(x.sp@bmat)>=1000]

bmat.t = bmat[x.sp@cluster==cluster& rowSums(x.sp@bmat)>=1000,]
bmat2 = bmat.t

#### subsample to 1000 reads. 
NUM=1000
for (i in 1:nrow(bmat2)) {
print(i)
if (length(which(bmat2[i,]==1)) > 1000) {
idx = sample(which(bmat2[i,]==1),length(which(bmat2[i,]==1))-NUM )
bmat2[i,idx] = 0
}
}

rsums = rowSums(bmat.t)

p1 = data.frame(samples,rsums)


calJaccard <- function(X_i, X_j){
  A = Matrix::tcrossprod(X_i, X_j);
  bi = Matrix::rowSums(X_i);
  bj = Matrix::rowSums(X_j);
  jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
  return(jmat);       
}


jcard = calJaccard(bmat2,bmat2)

rownames(jcard) = colnames(jcard) = samples

library(gdata)

jc.mat = unmatrix(jcard)

melted = data.frame(names=names(jc.mat),jc = jc.mat)

melted = melted[which(melted$jc!=1),]


#x = matrix(c(1,1,1,0,1,1,0,0,0,1,0,1),nrow=3)
#calJaccard(x,x)


agg = aggregate(jc~names,melted,mean)
#agg = aggregate(jc~names,melted,median)

agg$x = sub("(.*):(.*)","\\1",agg$names)
agg$y = sub("(.*):(.*)","\\2",agg$names)

system("mkdir ../epigenetic_drift")
pdf(paste0("../epigenetic_drift/",cluster,".epi_drift.pdf"))
ggplot(p1) + geom_boxplot(aes(samples,rsums))
ggplot(agg) + geom_tile(aes(x=x,y=y,fill=jc))
ggplot(melted) + geom_boxplot(aes(x=names,y=jc)) + 
  coord_flip()

dev.off()







