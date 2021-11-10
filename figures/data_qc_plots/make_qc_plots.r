#sample=commandArgs(trailing=T)[1]

sample="BM_03_rep1"
files = list.files(path="../../../data/snATAC/qc/bc_info_by_sample/")
samples = sub("(.*).barcode.info.txt","\\1",files)


dat_list = list()

for (sample in samples){ 

a=read.delim(paste0("../../../data/snATAC/qc/bc_info_by_sample/",sample,".barcode.info.txt"))
a$log10UM = log10(a$filter)
a = a[which(a$filter>=100),]
a$sample = sample
dat_list[[sample]] = a

}

dat = do.call(rbind,dat_list)
dat$tissue = substr(dat$sample,1,2)
dat$age_rep = substr(dat$sample,4,10)


pdf("allsample_tss_enrich.vs.counts.pdf",height=15,width=15)

ggplot(dat,aes(x=filter,y=TSS_enrich)) + 
  scale_x_log10() +
#  geom_point(size=0.1,alpha=0.1) + 
  geom_density_2d(alpha=0.5) +
  geom_density_2d_filled(alpha=0.5) +
  facet_wrap(tissue~age_rep)+
  ylim(0,25) + 
  scale_fill_manual(values=hsv(1, seq(0,1,length.out = 11) , 1)) +
  theme_bw()
dev.off()  

