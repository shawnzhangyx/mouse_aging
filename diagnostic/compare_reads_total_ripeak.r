tissue="heart"

b= read.table("../../../data/snATAC/ct2peaks/heart.all.barcode.peak.cnts")
b2 = read.table("../../../data/snATAC/bam_bowtie2_Olivier/qc/heart.all.dedup.filter.barcode.cnts")
b3 = read.table("../../../data/snATAC/bam_bowtie2_Olivier/qc/heart.all.barcode.cnts")
## cluster info
b4 = read.table("../../../analysis/Yang_NMF_method/heart/R6/heart.R6.statH")
dat = merge(b,b2,by="V1")
dat$ripr = dat$V2.x/dat$V2.y
dat$cluster = (b4$V3+1)[match(dat$V1,b4$V1)]
dat$cluster[is.na(dat$cluster)] = 0
ggplot(dat) + geom_point(aes(x=log10(V2.y),y=ripr,color=factor(cluster)),size=0.5)+
    facet_wrap(.~cluster)

                         