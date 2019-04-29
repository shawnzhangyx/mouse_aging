library(dplyr)
tissue="heart"
a = read.table("../../analysis/Yang_NMF_method/heart/heart.statH")

sample="03.rep1"

int.list = list()
for (sample in c("03.rep1","03.rep2","10.rep1","10.rep2","18.rep1","18.rep2")) {
print(sample)
b = fread(paste0("../../data/snATAC/bam_bowtie2_Olivier/insert_size/heart.",sample,".dedup.filter.insert_size.txt.gz"))
b$V1 = paste0(sample,".",b$V1)
b2 = b[which(b$V1 %in% a$V1),]
int.list[[sample]] = b2
}

dat = do.call(rbind,int.list)

dat2 = dat[which(dat$V2>0),]

dat2$cluster = a$V3[match(dat2$V1,a$V1)]
#dat3= subset(dat2, cluster==0)
dat3 = dat2 %>% count(V2,cluster)


pdf("cluster.vs.insert_size_dist.pdf",height=10,width=3)
ggplot(dat3) + geom_col(aes(x=V2,y=n)) + xlim(0,500) + 
  facet_grid(cluster~.,scales="free")
dev.off()

