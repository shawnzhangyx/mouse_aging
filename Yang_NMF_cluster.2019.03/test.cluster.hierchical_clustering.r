tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]

tissue="dorsal_hippocampus"
rank=15
setwd(paste0("../../analysis/Yang_NMF_method/",tissue,"/R",rank))
a=read.delim(paste0(tissue,".counts"),skip=1)
b =read.delim(paste0(tissue,".counts.summary"))
cnts = as.tibble(a[,-c(1:6)])
#total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])
grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..bam","\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)

cnts2 = t(cnts)
cnts3 = rowsum(cnts2,group=grps)
cnts3 = sweep(cnts3,1,rowSums(cnts3),"/")
cor.mat  = cor(t(cnts3))

hc = hclust(as.dist(1-cor.mat))
hc$order

labels_o = rev(c("C10","C15","C6","C2","C14","C5","C12","C3","C7","C9",
                   "C1","C8","C11","C13","C4"))

#hc$order = match(labels_o,hc$labels)
pdf(paste0(tissue,".R",rank,".hclust.pdf"))
plot(as.dendrogram(hc),horiz=T)
dev.off()





