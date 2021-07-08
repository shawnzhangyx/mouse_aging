setwd("../../analysis/paired_tag_h3k9me3/bam_split_cell")

b1 = read.delim("L23.03.rep1.L23_changed_domain.counts.summary")
b2 = read.delim("L23.03.rep1.read.counts.summary")
dat1 = as.data.frame(t(rbind(b1[1,-1],b2[1,-1],colSums(b1[,-1]))))
dat1$sample = "m3rep1"
dat1$age = "m3"

b1 = read.delim("L23.03.rep2.L23_changed_domain.counts.summary")
b2 = read.delim("L23.03.rep2.read.counts.summary")
dat2 = as.data.frame(t(rbind(b1[1,-1],b2[1,-1],colSums(b1[,-1]))))
dat2$sample = "m3rep2"
dat2$age = "m3"

b1 = read.delim("L23.18.rep1.L23_changed_domain.counts.summary")
b2 = read.delim("L23.18.rep1.read.counts.summary")
dat3 = as.data.frame(t(rbind(b1[1,-1],b2[1,-1],colSums(b1[,-1]))))
dat3$sample = "m18rep1"
dat3$age = "m18"

b1 = read.delim("L23.18.rep2.L23_changed_domain.counts.summary")
b2 = read.delim("L23.18.rep2.read.counts.summary")
dat4 = as.data.frame(t(rbind(b1[1,-1],b2[1,-1],colSums(b1[,-1]))))
dat4$sample = "m18rep2"
dat4$age = "m18"

dat = rbind(dat1,dat2,dat3,dat4)


colnames(dat)[1:3] = c("c.h3k9","h3k9","all")
dat$c.h3k9.over.h3k9 = dat$c.h3k9/dat$h3k9
dat$c.h3k9.over.all = dat$c.h3k9/dat$all
dat$h3k9.over.all = dat$h3k9/dat$all
dat$age = factor(dat$age,levels=c("m3","m18"))
#ggplot(dat) +geom_violin(aes(x=age,y=h3k9.over.all))

## change the maximum to 20%. 
dat$c.h3k9.over.all[dat$c.h3k9.over.all>0.2] = 0.2
#ggplot(dat) +geom_violin(aes(x=age,y=c.h3k9.over.all))

pdf("Fraction_of_reads_in_diff_h3k9me3_domain.FC.L23.pdf",width=4,height=4)
ggplot(subset(dat,all>200)) +geom_violin(aes(x=age,y=c.h3k9.over.all)) +
#  scale_x_discrete(breaks=c("m3","m18"))+ 
  theme_classic()
dev.off()

#ggplot(subset(dat,all>=200)) +geom_histogram(aes(x=c.h3k9.over.all),bins=100) +
#  scale_x_continuous(limits=c(0,0.1)) + 
#  facet_grid(age~.)



