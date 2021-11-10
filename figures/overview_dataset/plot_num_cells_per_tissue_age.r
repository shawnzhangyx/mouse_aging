
dat_list = list()
for (  tissue in c("DH","FC","HT","LM","BM") ){
 a =read.delim(paste0("../../../analysis/snapATAC/",tissue,"/",tissue,".pool.barcode.meta_info.txt"))
 dat_list[[tissue]] = a
}

dat = do.call(rbind,dat_list)

dat$tissue = substr(dat$sample,1,2)

meta=read.delim("../celltype_annotation.txt")
dat$Clade = meta$Clade[match(paste(dat$tissue,dat$cluster),paste(meta$Tissue,meta$Cluster))]
dat = dat[!is.na(dat$Clade),]

library(gridExtra)
g1=ggplot(dat) + geom_bar(aes(x=tissue,fill=factor(stage,levels=c(3,10,18))),stat="count",position="stack") + coord_flip() + scale_fill_brewer(palette="YlGnBu",name="Age(Months)") + 
  theme_classic()

g2=ggplot(dat) +geom_bar(aes(x=tissue,fill=factor(Clade)),stat="count",position="stack") + 
  coord_flip() + scale_fill_discrete("Clade") +
  theme_classic()

pdf("num_cells.pdf",width=8,height=6)
grid.arrange(g1,g2)

dev.off()

