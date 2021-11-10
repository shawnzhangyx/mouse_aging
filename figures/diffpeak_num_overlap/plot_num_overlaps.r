
a = read.table("overlaped.c.bed")
a = a[order(-a$V4),]
a$w = a$V3-a$V2
a = a[a$w<3000,]
a$tissue = NA
b = read.table("overlaped.wo.bed")
for (idx in 1:nrow(a)){
  print(idx)
  a$tissue[idx] = paste(b$V7[which(paste(b$V1,b$V2,b$V3) == paste(a$V1[idx],a$V2[idx],a$V3[idx]))],collapse=",")
  }

write.table(a, "overlaped.details.txt",row.names=F,sep="\t",quote=F)


id1 = which(a$V4>=5 & a$V4<10)
id2 = which(a$V4>=10)
a$V4[id1] = "5-9"
a$V4[id2] = ">=10"

shared_diff = data.frame(table(a$V4))
shared_diff$Var1 = factor(shared_diff$Var1,levels=c("1","2","3","4","5-9",">=10"))


pdf("All_tissues_Diff_peaks_shared_by_clusters.01.pdf",height=4,width=3)
ggplot(shared_diff) + geom_col(aes(x=Var1,y=Freq)) +
  geom_text(aes(x=Var1,y=Freq,label=Freq),nudge_y=1000)+
  xlab("# of Cell Types") + ylab("# of diff peaks") +
  theme_light()
dev.off()


