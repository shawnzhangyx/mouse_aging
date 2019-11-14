tissue="DH"
rank = 15

setwd(paste0("../../analysis/snapATAC/",tissue,"/snapFiles/"))
load(paste0(tissue,".pool.snapATAC.cluster.RData"))


sample.barcode = paste(x.sp@sample,x.sp@barcode,sep=":")

# read TSS enrichment scores. 

files = list.files(pattern="DH",path="/mnt/tscc/lamaral/projects/Aging/TSS_enrich",full.names=T)

tss_list = list()
for (file in files) {
  sample = sub(".*(DH_.._....)_.*","\\1",file)
  tss_list[[file]] =fread(file)
  tss_list[[file]]$sample = sample
}
tss_dat = do.call(rbind,tss_list)

tss_dat$sample.barcode = paste(tss_dat$sample,tss_dat$V1,sep=":")

tss_dat$cluster = x.sp@cluster[match(tss_dat$sample.barcode, sample.barcode)]

tss_dat2 = tss_dat[!is.na(tss_dat$cluster),]

pdf("../QC/TSS_enrichment_by_cluster.pdf")
ggplot(tss_dat2) + geom_violin(aes(cluster,V4)) +
geom_boxplot(aes(cluster,V4),width=0.2) + 
  geom_hline(yintercept=7,color="red",linetype="dashed")+
  xlab("cluster") + ylab("TSS enrichment")
dev.off()

ggplot(tss_dat2) + geom_density(aes(V4,fill=cluster==4),alpha=0.5) +
  geom_vline(xintercept=7,color="red",linetype="dashed")

