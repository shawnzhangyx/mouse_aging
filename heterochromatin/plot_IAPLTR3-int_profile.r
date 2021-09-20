setwd("../../analysis/heterochromatin")
a=read.table("FC.IAPLTR3-int.cluster_age.matrix.mat",skip=1)
b=readLines("FC.IAPLTR3-int.cluster_age.matrix.mat",n=1)

samples = unlist(strsplit(b,","))[393:446]
samples = sub(".*FC.metacell_(.*).sorted.rpkm.*","\\1",samples)
cluster = sub("(.*)\\.(.*)","\\1",samples)
age = sub("(.*)\\.(.*)","\\2",samples)

meta = read.delim("../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
meta = meta[which(meta$Tissue=="FC"),]
meta$Clade[9] = "Glia"
clade = meta$Clade[match(cluster,meta$Cluster)]

name = meta$Name[match(cluster,meta$Cluster)]

a2= a[,-c(1:6)]
mean = colSums(a2)/nrow(a2)

len = length(a2)/length(samples)

dat = data.frame(sample=rep(samples,each=len),
    pos = rep(1:len,length(samples)),
    mean =mean,
    age = rep(age,each=len),
    cluster= rep(cluster,each=len),
    clade = rep(clade,each=len),
    name = factor(rep(name,each=len))
    )

WIN =2 

dat2 = dat[which(dat$pos>WIN & dat$pos < 100-WIN),]

for ( sample in unique(dat2$sample)) {
  for ( pos in unique(dat2$pos)) {
    print(as.character(sample))
    print(pos)
    dat2$mean[which(dat2$sample==sample & dat2$pos==pos)] = mean(dat$mean[which(dat$sample==sample & dat$pos %in% (pos-WIN):(pos+WIN))],na.rm=T) 
    }
    }


pdf("FC.IAPLTR3.profile.pdf",height=8,width=16)
ggplot(dat2) + geom_line(aes(x=pos,y=mean,color=clade,group=sample)) +
  geom_point(data=subset(dat2,pos%%50==0), aes(x=pos,y=mean,shape=name,color=clade),size=5) + 
  scale_shape_manual(values=1:nlevels(dat2$name)) + 
  facet_wrap(.~age) + 
  theme_classic()
dev.off()


