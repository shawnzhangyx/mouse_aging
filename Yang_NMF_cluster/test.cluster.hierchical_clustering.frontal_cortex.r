tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]

tissue="frontal_cortex"
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

labels_o = rev(c("C12","C4","C2","C5","C6","C3","C7","C9","C11","C8",
                   "C13","C1","C10","C14","C15"))

hc$order = match(labels_o,hc$labels)
pdf(paste0(tissue,".R",rank,".hclust.pdf"))
plot(hc)
#plot(as.dendrogram(hc),horiz=T)
dev.off()





