
a = read.table("overlaped.clades.bed")
a[,4:10][a[,4:10]>1] = 1
a$total = rowSums(a[,4:10])
#a = a[order(-a$V4),]
a$w = a$V3-a$V2
a = a[a$w<3000,]
a = a[a$total>0,]

shared_diff = data.frame(table(a$total))

pdf("All_tissues_Diff_peaks_shared_by_clades.01.pdf",height=4,width=3)
ggplot(shared_diff) + geom_col(aes(x=Var1,y=Freq)) +
  geom_text(aes(x=Var1,y=Freq,label=Freq),nudge_y=1000)+
  xlab("# of Cell Classes") + ylab("# of diff peaks") +
  theme_light()
dev.off()


