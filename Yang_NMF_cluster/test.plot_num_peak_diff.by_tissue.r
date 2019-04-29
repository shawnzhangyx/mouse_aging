
summ1 = read.delim("../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/age_diff_bycluster/dorsal_hippocampus.peaks.age_diff.summary.txt")
summ2 = read.delim("../../analysis/Yang_NMF_method/heart/R10/age_diff_bycluster/heart.peaks.age_diff.summary.txt")
summ3 = read.delim("../../analysis/Yang_NMF_method/frontal_cortex/R15/age_diff_bycluster/frontal_cortex.peaks.age_diff.summary.txt")

out = data.frame(DH=colSums(summ1[,-1]),HT=colSums(summ2[,-1]),FC=colSums(summ3[,-1]))

melted = melt(as.matrix(out))
melted$Var2 = factor(melted$Var2, levels=c("DH","FC","HT"))

pdf("../../analysis/Yang_NMF_method/all_tissue.age_diff.summary.pdf",width=5,height=3)
ggplot(melted) + geom_col(aes(x=Var2,y=value,fill=Var1),position="stack") + 
  coord_flip() + theme_bw() 
dev.off()                                                  