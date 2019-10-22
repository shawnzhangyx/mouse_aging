library(Matrix)

## read the matrix of repeats
# row is cell, column is feature(repeat).
a=readMM("DH_03_rep1.mtx")
x=readLines("DH_03_rep1.xgi")
y=readLines("DH_03_rep1.ygi")
a1 = as.matrix(a)
#rownames(a1) = x
#colnames(a1) = y
## read the total number of reads.
b=read.delim("../../data/snATAC/bam.filter/DH_03_rep1/DH_03_rep1.filter.barcode.cnts.txt",header=F)
## read the cluster information
d=read.delim("../../analysis/snapATAC/DH/DH.pooled.barcode.cluster.stage.rep.txt")
d = d[which(d$stage==3 & d$rep=="rep1"),]

a2 =sweep(a1,1,b$V2[match(x,b$V1)],'/')
rownames(a2) = x
colnames(a2) = y
a3 = a2[which(x %in% d$barcode),]
#cluster = d$cluster[match(rownames(a3),d$barcode)]

melted = melt(a3)
melted$cluster = d$cluster[match(melted$Var1,d$barcode)]

avg = aggregate(value~Var2+cluster,melted,mean)


