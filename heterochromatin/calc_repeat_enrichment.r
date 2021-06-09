setwd()
a=read.table("overlap_repeats.counts.txt",header=F)
a = a[,c(1,2,4,6,7)]

a= a[order(-a$V6),]

a$pval = apply(a[,c(2:5)], 1,function(vec){ fisher.test(matrix(vec,nrow=2))$p.value })
a$fdr = p.adjust(a$pval,method="BH")
a = a[order(a$fdr),]
a$fc = (a$V6/a$V7)/(a$V2/a$V4)


a$V1 = factor(a$V1,levels=rev(a$V1))

pdf("enriched_repeats.pdf",height=2,width=3)
ggplot(head(a)) + geom_col(aes(x=V1,y=-log10(fdr))) +
  coord_flip()
dev.off()



