setwd("../../analysis/snapATAC/all_celltypes")
# read the order of cell types
od = read.delim("Spearman.clustering.celltypes.ordered.txt",
    header=F,stringsAsFactors=F)$V1
ct = sub("(..)\\.(.*?)\\..*","\\1.\\2",od)

## num of reads per cluster
reads= read.delim("all_celltypes.counts.summary")
reads = data.frame(colSums(reads[,-1]))
colnames(reads) = "reads"
reads$ct = sub(".*\\.(..)\\.metacell_(.*)\\.sorted.bam","\\1.\\2",rownames(reads))

reads = reads[match(ct,reads$ct),]

## num of cells per cluster
dat_list = list()
for (  tissue in c("DH","FC","HT","LM","BM") ){
 a =read.delim(paste0("../../../analysis/snapATAC/",tissue,"/",tissue,".pool.barcode.meta_info.txt"))
  dat_list[[tissue]] = a
  }
dat = do.call(rbind,dat_list)
dat$tissue = substr(dat$sample,1,2)
dat$ct = paste0(dat$tissue,".",dat$cluster)

cellnum = data.frame(table(dat$ct))

cellnum = cellnum[match(ct,cellnum$Var1),]

# combine the columns. 
a = data.frame(od=factor(od,levels=od),reads=reads$reads,cellnum = cellnum$Freq)


g1 = 
  ggplot(a) + geom_col(aes(x=od,y=cellnum)) +
  geom_text(aes(x=od,y=cellnum,label=cellnum),hjust=-0.2) + 
  coord_flip() + ylim(0,15000) +
  theme_bw()

g2 = 
  ggplot(a) + geom_col(aes(x=od,y=reads/1e6)) + 
  geom_text(aes(x=od,y=reads/1e6,label=paste0(format(reads/1e6,digits=2),"M")),hjust=-0.2) +
  coord_flip() + scale_y_continuous(name="reads(million)",limits=c(0,1e2))+
  theme_bw()# + 
#  theme(axis.text.y = element_blank())


library(gridExtra)
pdf("all_celltypes.cell_num_reads.pdf",width=10,height=15)
grid.arrange(g1,g2,nrow=1)
dev.off()




