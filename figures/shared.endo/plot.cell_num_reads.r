setwd("../../../analysis/snapATAC/all_celltypes")
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
#dat$ct.stage.rep = paste0(dat$ct,".",dat$stage,".",dat$rep)



cellnum = aggregate(sample~ct+stage+rep,dat,length) #data.frame(table(dat$ct))

cellnum$od = od[match(cellnum$ct,ct)]

a= cellnum[grepl("Endo",cellnum$od),]

b = aggregate(sample~od,a,sum)
b = b[order(b$sample),]
b$od = factor(b$od,levels=b$od)
#a=a[order(a$stage),]

g1 = 
  ggplot(a) + geom_col(aes(x=od,y=sample,fill=paste0(stage,".",rep)),position="dodge") +
#  geom_text(aes(x=od,y=cellnum,label=cellnum),hjust=-0.2) + 
#  coord_flip() + ylim(0,15000) +
  theme_bw()

g2 = 
  ggplot(b) +  geom_col(aes(x=1,y=sample,fill=od),position="stack") +
  theme_classic()


setwd("../../../aging_share/figures/shared.endo/")

pdf("endo.cell_num.tissue.pdf",width=3,height=5)
g2
dev.off()

pdf("endo.cell_num.tissue_age_rep.pdf",width=10,height=5)
g1
dev.off()



