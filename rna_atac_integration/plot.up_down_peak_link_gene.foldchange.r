
setwd("../../analysis/rna_atac_integration")

up.dh = read.csv("DH.up.out.csv")
up.fc = read.csv("FC.up.out.csv")
up_out = rbind(up.dh,up.fc)

down.dh = read.csv("DH.down.out.csv")
down.fc = read.csv("FC.down.out.csv")
down_out = rbind(down.dh,down.fc)




pdf("logFC_up_down_peaks.combined.pdf",width=4,height=2)
ggplot(up_out) + #geom_boxplot(aes(x=1,y=logFC),outlier.shape=NA) +
  geom_boxplot(data=up_out, aes(x="up",y=logFC),outlier.shape=NA,width=0.3) +
  geom_boxplot(data=down_out,aes(x="down",y=logFC),outlier.shape=NA,width=0.3) +
  geom_hline(yintercept=0,linetype="dashed") +
  ylim(-5,5) + 
  coord_flip() + theme_bw()
dev.off()

w = wilcox.test(up_out$logFC,down_out$logFC)
# 1.8e-13
