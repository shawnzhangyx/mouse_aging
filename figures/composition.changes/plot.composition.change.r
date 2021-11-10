tissue="DH"

ct = read.delim("../celltype_annotation.txt")

for (tissue in c("DH","FC","HT","LM","BM")) {
  print(tissue)
  a=read.delim(paste0("../../../analysis/snapATAC/",tissue,"/",tissue,".pool.barcode.meta_info.txt"))
  a$stage[a$stage==3] = "03"
  a$sample = paste(a$stage,a$rep)
  a$ct = paste0(ct$Tissue,".",ct$Cluster,".",ct$Name)[match( paste(tissue,a$cluster), paste(ct$Tissue,ct$Cluster))]
  a$ct = factor(a$ct, levels=paste0(ct$Tissue,".",ct$Cluster,".",ct$Name)[which(ct$Tissue==tissue)])

tab = table(a$ct,a$sample)
tab = sweep(tab,2,colSums(tab),'/')


melted = melt(tab)
g2 = ggplot(melted) +
  geom_col(aes(factor(Var1),value,
  fill = factor(Var2)
  ),color="black",
  position="dodge") +
  scale_fill_manual(values=c("grey1","grey2","grey31","grey32","grey61","grey62")) + 
  xlab("Cell Types") + ylab("Fraction") +
  guides(fill=guide_legend(title="Age Group"))  +
  theme_bw(base_size=10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1) )



pdf(paste0(tissue,".celltype.fraction.pdf"),height=4,width=8)
  print( g2)
  dev.off()

}


