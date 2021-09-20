setwd("../../analysis/all_celltypes.all_tissue")

a=read.delim("all_celltypes.thresh.counts",skip=1)

b= a[,-c(1:6)]
rownames(b) = a[,1]
b = sweep(b,2,colSums(b),'/') *1e6

b2 = b[which(rowSums(b)!=0),]

write.csv(b2,"all_celltypes.thresh.cpms.csv")
