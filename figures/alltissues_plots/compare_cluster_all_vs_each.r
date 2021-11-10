a=read.delim("luisa_all_tissues/all_merged_40kL.cluster.meta.txt")
m1 = read.csv("cluster_keys.csv")
m2 = read.delim("../celltype_annotation.txt")
m2 = m2[!is.na(m2$Clade),]

## correct LM cells.
LM = read.delim("../../../analysis/snapATAC/LM/LM.pool.barcode.meta_info.txt")
LM$tissue_cluster = paste0("LM_",LM$cluster)
tmp = a[which(a$tissue=="LM"),]
a[which(a$tissue=="LM"),"tissue_cluster"] = LM$tissue_cluster[match(paste(tmp$sample,tmp$barcode),paste(LM$sample,LM$barcode))]


a$ct_all = m1$ct[match(a$cluster,m1$cluster)]
a$ct_tissue = paste0(m2$Tissue,".",m2$Name)[match(a$tissue_cluster,paste0(m2$Tissue,"_",m2$Cluster))]

tab = table(a$ct_all,a$ct_tissue)
#write.csv(tab,"all.vs.tissue.cluster.csv")

#sort the table. 
tab2 = as.data.frame.matrix(tab)
tab2 = tab2[order(-rowSums(tab2)),]

j=1
for (i in 1:nrow(tab2)) {
  #od = order(!ct.max==i)
  od = order(-tab2[i,j:ncol(tab2)])
  tab2[,j:ncol(tab2)] = tab2[,j:ncol(tab2)][,od]
  colnames(tab2)[j:ncol(tab2)] = colnames(tab2)[j:ncol(tab2)][od]
  ct.max = apply(tab2[j:ncol(tab2)],2,which.max)
  j = j + length(which(ct.max==i))
}

tab2[1:5,1:5]

write.csv(tab2,"all.vs.tissue.cluster.csv")

tab3 = sweep(tab2,2, colSums(tab2),'/')
library(pheatmap)
pdf("all.vs.tissue.cluster.heatmap.pdf",height=10,width=30)
pheatmap(tab2, cluster_rows =F, cluster_cols=F,
    color=colorRampPalette(c("white","red"))(50),
    display_numbers=T,number_format="%i",fontsize=15)
dev.off()
pdf("all.vs.tissue.cluster.heatmap.pc.pdf",height=10,width=30)
pheatmap(tab3, cluster_rows =F, cluster_cols=F,
  color=colorRampPalette(c("white","red"))(50),
  display_numbers=T,fontsize=15)
dev.off()


tab4 = sweep(tab2,1, rowSums(tab2),'/')
pdf("all.vs.tissue.cluster.heatmap.pc_row.pdf",height=10,width=30)
pheatmap(tab4, cluster_rows =F, cluster_cols=F,
  color=colorRampPalette(c("white","red"))(50),
  display_numbers=T,fontsize=15)
dev.off()

sumP = 0 

#j=1
for (i in 1:nrow(tab2)) {
  ct.max = apply(tab2,2,which.max)
  rowSum = sum(tab2[i,])
  if (length(which(ct.max==i)) >0) {
  ctSum = sum(tab2[i,which(ct.max==i)]) }
  else { ctSum=0 } 
  sumP[i] = ctSum/rowSum
}

write.csv(data.frame(rownames(tab2),sumP),"major_celltype_fraction.csv")




