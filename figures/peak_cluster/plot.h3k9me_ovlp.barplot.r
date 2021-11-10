a=read.table("FB.P0.H3K9me3.gt100k_domain.ovlp_ATACpeaks.bed")
#b = a[which(a$V8>0),]
b=read.table("clustered_diff_peaks.merged.bed")

GS = 2730871774

up = sum((b$V3-b$V2)[which(b$V4=="Up")])
down = sum((b$V3-b$V2)[which(b$V4=="Down")])


#num = c(sum(a$V5)/GS,sum(a$V10)/GS, sum(a$V5)/GS*up/GS,0,sum(a$V5)/GS*down/GS)
#names = c("H3K9me3","ovlp.Up.Obs","Ovlp.Up.Exp","Ovlp.Down.Obs","Ovlp.Down.Exp")
num = c(sum(a$V5)/GS,sum(a$V10)/GS, sum(a$V5)/GS*up/GS)*100
names = c("H3K9me3","ovlp.Up.Obs","Ovlp.Up.Exp")

dat = data.frame(names=factor(names,levels=names),num=num)

pdf("h3k9me3.ovlp_w_peak_cluster.pdf",height=3,width=2.5)
ggplot(dat) + geom_bar(aes(x=names,y=num),stat="identity",width=0.6) +
  geom_text(aes(x=names,y=num,label=format(num,digits=2)),vjust=-1) + 
  scale_y_continuous(limits=c(0,3.5)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1) 
  )
dev.off()
