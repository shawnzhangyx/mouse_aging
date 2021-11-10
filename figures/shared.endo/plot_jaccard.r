library(pheatmap)
library(viridis)
a=read.csv("../diffpeak_num_overlap/All_tissues.diffpeak_share.table.1pc.csv",row.names=1)

mat = a[grepl("Endo",rownames(a)),grepl("Endo",colnames(a))]

MAX  = max(mat[mat!=1])*2# 0.02

mat[which(mat==1, arr.ind = T)] = MAX
#mat[which(mat>MAX, arr.ind = T)] = MAX


pdf("endo.jaccard.1pc.pdf",height=8,width=9)
pheatmap(mat,color=viridis(50),fontsize=20)
dev.off()

a=read.csv("../diffpeak_num_overlap/All_tissues.diffpeak_share.table.p0.01.csv",row.names=1)
mat = a[grepl("Endo",rownames(a)),grepl("Endo",colnames(a))]
MAX = max(mat[mat!=1])*2 # 0.02
mat[which(mat==1, arr.ind = T)] = MAX
mat[which(mat>MAX, arr.ind = T)] = MAX


pdf("endo.jaccard.p0.01.pdf")
pheatmap(mat)
dev.off()




a=read.csv("../diffpeak_num_overlap/All_tissues.diffpeak_share.table.top1k.csv",row.names=1)
mat = a[grepl("Endo",rownames(a)),grepl("Endo",colnames(a))]
MAX = max(mat[mat!=1])*2#0.02
mat[which(mat==1, arr.ind = T)] = MAX
mat[which(mat>MAX, arr.ind = T)] = MAX


pdf("endo.jaccard.p0.01.pdf")
pheatmap(mat)
dev.off()



