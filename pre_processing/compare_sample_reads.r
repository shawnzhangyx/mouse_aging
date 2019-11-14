setwd("../../data/snATAC/")

cnts = read.delim("counts/DH.read.counts",skip=1)
summ = read.delim("counts/DH.read.counts.summary")
sums = colSums(summ[,-1])

colnames(cnts)[-c(1:6)] = sub(".*.(.._.._rep.).*","\\1",colnames(cnts)[-c(1:6)])

cnts2 = cnts[-c(1:6)] 
rownames(cnts2) = cnts$Geneid
rpm = sweep(cnts2,2,colSums(cnts2),"/")

arpm = sweep(cnts2,2,sums,"/")


#ggplot(rpm) + geom_point(aes(DH_03_rep1,DH_03_rep2),color="grey") +
#  geom_abline(slope=1) + 
#  scale_x_log10() + scale_y_log10()

#ggplot(rpm) + geom_point(aes(DH_03_rep1,DH_18_rep1),color="grey") +
#  geom_abline(slope=1) +
#  scale_x_log10() + scale_y_log10()


pdf("qc/compare_DH_samples_with_different_TSS_enrich_scores.pdf")
ggplot(rpm) + geom_bin2d(aes(log10(DH_03_rep1+1e-5),log10(DH_03_rep2+1e-5)),bins=200)+
  geom_abline(slope=1,color='red') + ggtitle("Normalize by Reads in Peaks")

ggplot(arpm) + geom_bin2d(aes(log10(DH_10_rep1+1e-5),log10(DH_10_rep2+1e-5)),bins=200)+
  geom_abline(slope=1,color='red') + ggtitle("Normalize by Total Reads")

ggplot(rpm) + geom_bin2d(aes(log10(DH_10_rep1+1e-5),log10(DH_10_rep2+1e-5)),bins=200) +
  geom_abline(slope=1,color='red') + ggtitle("Normalize by Reads in Peaks")


ggplot(arpm) + geom_bin2d(aes(log10(DH_10_rep1+1e-5),log10(DH_10_rep2+1e-5)),bins=200)+
  geom_abline(slope=1,color='red') + ggtitle("Normalize by Total Reads")


dev.off()

