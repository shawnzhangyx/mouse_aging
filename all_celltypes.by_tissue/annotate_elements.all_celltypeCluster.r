setwd("../../analysis/snapATAC/all_celltypes")
#H3K27ac H3K27me3 H3K4me3 H3K4me1 H3K9me3 CTCF

a=data.frame(fread("allcelltypes.cluster_peaks.merged.annotated.txt"))

a2  = a[,4:9]

colnames(a2) = c("H3K27ac","H3K27me3","H3K4me3","H3K4me1","H3K9me3","CTCF")


state = rep("NS",nrow(a2))
state[which(a2$H3K4me1>0)] = "Enh-InA"
state[which(a2$H3K4me3>0)] = "Pr-InA"

state[which(a2$H3K27ac>0 & a2$H3K4me1>0)] = "Enh-A"
state[which(a2$H3K27ac>0 & a2$H3K4me3>0)] = "Pr-A"

state[which(a2$H3K27me3>0)] = "K27me-Hetero"
state[which(a2$H3K27me3>0 & a2$H3K4me1>0)] = "Enh-Bi"
state[which(a2$H3K27me3>0 & a2$H3K4me3>0)] = "Pr-Bi"

state[which(a2$CTCF>0)] = "CTCF"
state[which(a2$H3K9me3>0)] = "K9-Hetero"


tab = as.data.frame(table(state))
tab$Ratio = tab$Freq/sum(tab$Freq)
tab$state = factor(tab$state, levels=rev(c("NS","CTCF","Enh-A","Enh-InA","Enh-Bi","Pr-A","Pr-InA","Pr-Bi","K27me-Hetero","K9-Hetero")))

pdf("annotation_of_peaks.pdf",height=4,width=4)
ggplot(tab) + geom_bar(aes(x=state,y=Ratio*100,fill=state),stat="identity",width=0.8) +
#    geom_hline(yintercept=0,size=0.5) + 
    scale_y_continuous(
    name = "%Peaks",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*sum(tab$Freq)/100, name="#Peaks")
  ) +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position="none",
#  panel.grid = element_blank() 
  ) 
dev.off()


