
library(gridExtra)

files = list.files(pattern="barcode.info",path="/projects/ps-renlab/yanxiao/projects/mouse_aging/data/snATAC/qc/bc_info_by_sample",full.names=T)

dat = list()
for (file in files){
tmp = fread(file)
tmp$sample = sub(".*\\/(.._.._rep.).barcode.info.txt","\\1",file)
dat[[file]] = tmp
}

dat2 = data.frame(do.call(rbind,dat))
dat2$tissue = substr(dat2$sample,1,2)
dat2$stage = substr(dat2$sample,4,5)
dat2$rep = substr(dat2$sample,7,10)



g1 = ggplot(subset(dat2,filter>=500)) +
  geom_bar(aes(sample,fill=tissue,alpha=TSS_enrich>=7),color="gray",stat="count",position="stack")+
  geom_hline(yintercept=5000,color="red",linetype="dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g2 = ggplot(subset(dat2,filter>=500)) +
  geom_bar(aes(sample,fill = TSS_enrich>=7),stat="count",position="fill")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

g3 = ggplot(subset(dat2,filter>=500)) +
  geom_bar(aes(sample,fill=tissue,alpha=TSS_enrich>=10),color="gray",stat="count",position="stack")+
  geom_hline(yintercept=5000,color="red",linetype="dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g4 = ggplot(subset(dat2,filter>=500)) +
  geom_bar(aes(sample,fill = TSS_enrich>=10),stat="count",position="fill")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#pdf("Number_Of_Usable_cells_Fragment_TSS_cutoff.pdf",height=10,width=8)
#grid.arrange(g1,g2,nrow=2)
#grid.arrange(g3,g4,nrow=2)
#dev.off()

g1 = ggplot(subset(dat2,filter>=100)) +
  geom_violin(aes(sample,TSS_enrich,fill=tissue),alpha=0.5) +
  geom_boxplot(aes(sample,TSS_enrich,fill=tissue),width=0.1,alpha=0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


g2 = ggplot(subset(dat2,filter>=100)) +
  geom_violin(aes(sample,log10(filter),fill=tissue),alpha=0.5) +
  geom_boxplot(aes(sample,log10(filter),fill=tissue),width=0.1,alpha=0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("Distrubtion_of_Fragment_TSS.pdf",height=10,width=10)
grid.arrange(g1,g2,nrow=2)
dev.off()



