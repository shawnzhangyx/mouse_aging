setwd('../../data/snATAC/bam_bowtie2_Olivier/qc/')
#tissue="dorsal_hippocampus"
tissue= commandArgs(trailing=T)[1]
names=list.files(pattern=paste0(tissue,".*cnts$"))
#names1 = names[!grepl("dedup",names)]
#names2 = names[grepl("dedup",names)]

tmp_list = {}
for (name in names){
a=read.table(name)
a$age = sub(".*.(..)\\.(....)(.picard.dedup)?\\.barcode.cnts","\\1",name)
a$rep = sub(".*.(..)\\.(....)(.picard.dedup)?\\.barcode.cnts","\\2",name)
a$dedup = ifelse(sub(".*.(..)\\.(....)(.picard.dedup)?\\.barcode.cnts","\\3",name)==".picard.dedup","yes","no")
a$sample = sub(".*.(..)\\.(....)(.picard.dedup)?\\.barcode.cnts","\\1.\\2.\\3",name)
tmp_list[[name]] = a
}

dat = do.call(rbind,tmp_list)

pdf(paste0(tissue,".barcode.cnts.pdf"))
ggplot(subset(dat,dedup=="no")) + geom_histogram(aes(log10(V2)),bins=50) + 
  facet_wrap(rep~age) + ggtitle("before dedup") + 
  theme_bw()
ggplot(subset(dat,dedup=="yes")) + geom_histogram(aes(log10(V2)),bins=50) +
  facet_wrap(rep~age) + ggtitle("after dedup") + 
  theme_bw()

ggplot(dat) + geom_violin(aes(x=sample,y=log10(V2),fill=rep,color=dedup)) +
  theme_bw()
dev.off()

#agg = aggregate(dat, V2~sample, function(x){length(which(x >2000)) 

agg = aggregate(V2~sample,dat, sum)
agg[seq(2,12,2),2]/agg[seq(1,12,2),2]
#[1] 0.8088344 0.5000476 0.7963848 0.4851282 0.7935821 0.4944541

