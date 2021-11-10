a=read.delim("clustered_diff_peaks.clustered.ovlpH3K9me3.bed",header=F)
a=a[which(a$V1!="chrY"),]

# column 5 peak cluster. 
# column 6 overlap with H3K9me3

a$ct = sub("(.*).(Up|Down)","\\1",a$V4) 
a$dir = sub("(.*).(Up|Down)","\\2",a$V4)

meta = read.delim("../celltype_annotation.txt")
meta$ct =paste0(meta$Tissue,".",meta$Cluster) 

a$clade = meta$Clade[match(a$ct,meta$ct)]

a2 = a[!duplicated(a$V5),]
a2$ovlp = a2$V6>0

table(a2$dir,a2$V6>0)

pdf("peak_cluster_ovlp_w_h3k9me3.pdf",height=3,width=3)
ggplot(a2) + geom_bar(aes(x=dir,fill=ovlp),position="stack",color="black",width=0.6) + 
  scale_fill_manual(values=c("grey","black")) +
  theme_classic()
dev.off()


