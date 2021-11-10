#tab = read.csv(paste("OR_tables/", "H3K9me3","_OR_tab.txt" , sep = ""), header = T, sep = "\t")

tab_list = list()
for (mark in c( "H3K27ac","H3K4me1","H3K27me3","H3K9me3","CTCF","H3K4me3")) {
  tab_list[[mark]] = read.csv(paste0("OR_tables/",mark,"_OR_tab.txt"),header = T, sep = "\t")
  tab_list[[mark]]$file = NULL
  tab_list[[mark]]$mark = mark
}

tab = do.call(rbind,tab_list)
meta = read.delim("../../celltype_annotation.txt")
tab$Name = meta$Name[match(paste(tab$Tissue,tab$Cluster),paste(meta$Tissue,meta$Cluster))]

tab$Clade = meta$Clade[match(paste(tab$Tissue,tab$Cluster),paste(meta$Tissue,meta$Cluster))]

#dd = dd[-which(dd$id %in% c("HT.doublet","HT.Low_quality","LM.Low_quality","FC.Inh.Lrrc38.Plcxd3.Prok2","LM.Eryth")),]
tab = tab[-which(paste0(tab$Tissue,".",tab$Name) %in% c("HT.doublet","HT.Low_quality","LM.Low_quality","FC.InN.Lrrc38.Plcxd3","LM.Eryth","DH.Low_quality")),]

tab$p[which(tab$OR<1)] = 1

tab$log10P = -log10(tab$p)
tab$mark = factor(tab$mark, levels=c("CTCF","H3K27ac","H3K4me3","H3K4me1","H3K27me3","H3K9me3"))
tab$Clade = factor(tab$Clade,levels=c("ExN","InN","Glia","Muscle","Myeloid","Lymphoid","Other"))
tab$id = paste0(tab$Tissue,".",tab$Name)
# sort each cell type by their sum of significance. 
agg.up = aggregate(log10P~id,subset(tab,dir=="Up"),sum)
agg.down = aggregate(log10P~id,subset(tab,dir=="Down"),sum)
agg = merge(agg.up,agg.down,by="id")
agg$diff = agg$log10P.x-agg$log10P.y
agg = agg[order(agg$diff),]
tab$id = factor(tab$id, levels=agg$id)

tab2 = tab 
MAX=20
tab2$log10P[which(tab2$log10P>MAX)] = MAX
#tab2$OR[which(tab2$OR>6)] = 6

pdf("diff_peak.histone_enrichment.pdf",height=16,width=6)

ggplot() + 
  geom_point(data=subset(tab2,log10P>0),aes(x=mark,y=id,size=log10P,fill=OR),shape=21) + 
  geom_text(data=subset(tab2,log10P==0),aes(x=mark,y=id,label="."))+ 
#  geom_point(aes(x=mark,y=paste(Tissue,Name,sep="."),fill=log10P,size=OR),shape=21) +
  scale_fill_gradient2(high="red",mid="white",low="white") + 
#  scale_fill_gradientn(colors=c("white","white","orange","red","red")) + 
  facet_grid(Clade~dir,scales="free_y",space="free_y") +
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1)) 

dev.off()


