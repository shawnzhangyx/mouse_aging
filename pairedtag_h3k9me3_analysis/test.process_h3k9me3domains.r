#a=read.table("peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland")
#a$w = a$V3-a$V2
#a = a[which(a$V1!="chrY"),]

a=read.delim("SICER-W5000-G10000.read.counts",skip=1)
colnames(a)[-c(1:6)] = sub("bam_merge_rep.(.*).bam","\\1",colnames(a)[-c(1:6)])

a[,-(1:6)] = sweep(a[,-(1:6)],2,colSums(a[,-(1:6)]),"/") *1e6

a$w = a$End-a$Start
a2 = subset(a,w>=1e5)


g1 = ggplot(a2)+ geom_point(aes(x=ASC.03,y=ASC.18),color="grey") + 
  geom_abline(intercept=0,slope=1) + 
  scale_x_log10()+
  scale_y_log10()


g2 = ggplot(a2)+ geom_point(aes(x=L23.03,y=L23.18),color="grey") +
  geom_abline(intercept=0,slope=1) +
  scale_x_log10()+
  scale_y_log10()

g3 = ggplot(a2)+ geom_point(aes(x=L4.03,y=L4.18),color="grey") +
  geom_abline(intercept=0,slope=1) +
  scale_x_log10()+
  scale_y_log10()


#library(gridExtra)
grid.arrange(g1,g2,g3)


g1 = ggplot(a2)+ geom_point(aes(x=log2(ASC.03),y=log2(ASC.18/ASC.03)),color="grey") +
  geom_hline(yintercept=0)

g2 = ggplot(a2)+ geom_point(aes(x=log2(L23.03),y=log2(L23.18/L23.03)),color="grey") +
  geom_hline(yintercept=0)

g3 = ggplot(a2)+ geom_point(aes(x=log2(L4.03),y=log2(L4.18/L4.03)),color="grey") +
  geom_hline(yintercept=0)


