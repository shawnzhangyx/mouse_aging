library(RColorBrewer)
library(gridExtra)
#library(ggsignif)
library(dplyr)

tissue = commandArgs(trailing=T)[1] 

setwd(paste0("/mnt/tscc/lamaral/projects/Aging/",tissue,"/overlap/"))
p = list()
clusters = length(list.files(path="peaks.bed",pattern="*.peaks.bed"))

for (i in 1:clusters) {
  H3K4me3 = read.csv(paste("peaks.bed/",i,".H3K4me3.loj.bed", sep = ""), sep = "\t", header = F) 
  head(H3K4me3)
  H3K4me3$id = paste0(H3K4me3$V1,":",H3K4me3$V2,"-",H3K4me3$V3, sep = "")
  if (length(which(duplicated(H3K4me3$id)))>0) {
    H3K4me3 = H3K4me3[-which(duplicated(H3K4me3$id)), ]
  }
  
  H3K4me1 = read.csv(paste("peaks.bed/",i,".H3K4me1.loj.bed", sep = ""), sep = "\t", header = F) 
  head(H3K4me1)
  H3K4me1$id = paste0(H3K4me1$V1,":",H3K4me1$V2,"-",H3K4me1$V3, sep = "")
  if (length(which(duplicated(H3K4me1$id)))>0) {
    H3K4me1 = H3K4me1[-which(duplicated(H3K4me1$id)), ]
  }
  
  
  CTCF = read.csv(paste("peaks.bed/",i,".CTCF.loj.bed", sep = ""), sep = "\t", header = F) 
  head(CTCF)
  CTCF$id = paste0(CTCF$V1,":",CTCF$V2,"-",CTCF$V3, sep = "")
  if (length(which(duplicated(CTCF$id)))>0) {
    CTCF = CTCF[-which(duplicated(CTCF$id)), ]
  }
  
  H3K27me3 = read.csv(paste("peaks.bed/",i,".H3K27me3.loj.bed", sep = ""), sep = "\t", header = F) 
  head(H3K27me3)
  H3K27me3$id = paste0(H3K27me3$V1,":",H3K27me3$V2,"-",H3K27me3$V3, sep = "")
  if (length(which(duplicated(H3K27me3$id)))>0) {
    H3K27me3 = H3K27me3[-which(duplicated(H3K27me3$id)), ]
  }
  
  
  H3K9me3 = read.csv(paste("peaks.bed/",i,".H3K9me3_reprocessed.loj.bed", sep = ""), sep = "\t", header = F) 
  head(H3K9me3)
  H3K9me3$id = paste0(H3K9me3$V1,":",H3K9me3$V2,"-",H3K9me3$V3, sep = "")
  if (length(which(duplicated(H3K9me3$id)))>0) {
    H3K9me3 = H3K9me3[-which(duplicated(H3K9me3$id)), ]
  }
  
  H3K27ac = read.csv(paste("peaks.bed/",i,".H3K27ac.loj.bed", sep = ""), sep = "\t", header = F) 
  head(H3K27ac)
  H3K27ac$id = paste0(H3K27ac$V1,":",H3K27ac$V2,"-",H3K27ac$V3, sep = "")
  if (length(which(duplicated(H3K27ac$id)))>0) {
    H3K27ac = H3K27ac[-which(duplicated(H3K27ac$id)), ]
  }
  
  all = data.frame(cbind(paste(H3K4me3$id), paste(H3K4me3$V7),paste(H3K27ac$V7),paste(H3K27me3$V7),paste(CTCF$V7),paste(H3K4me1$V7)), paste(H3K9me3$V7))
  colnames(all) = c("id","H3K4me3", "H3K27ac", "H3K27me3", "CTCF","H3K4me1", "H3K9me3")
  head(all)
  all$res = "NA"
  all$res[which(all$H3K4me3 != ".")] = "Prom"
  all$res[which(all$H3K4me3 == "." & all$H3K4me1 != ".")] = "Enh"
  all$res[which(all$H3K4me3 == "." & all$H3K4me1 == "." & all$CTCF != ".")] = "CTCF"
  all$res[which(all$H3K4me3 == "." & all$H3K4me1 == "." & all$CTCF == "." & all$H3K27me3 != ".")] = "Het"
  all$res[which(all$res == "Prom" & all$H3K27ac != ".")] = "Active Prom"
  all$res[which(all$res == "Enh" & all$H3K27ac != ".")] = "Active Enh"
  all$res[which(all$H3K4me3 == "." & all$H3K4me1 == "." & all$CTCF == "." & all$H3K27me3 == "." & all$H3K9me3 != ".")] = "Het"
  
  #pie(table(all$res))
  
  all$sig = paste0("NS-",nrow(all))
  pval_cut = sort(H3K4me3$V5)[nrow(H3K4me3)*0.01]
  all$sig[which(H3K4me3$V5<pval_cut & H3K4me3$V4<0)] = paste0("DOWN-",length(which(H3K4me3$V5<pval_cut & H3K4me3$V4<0)))
  all$sig[which(H3K4me3$V5<pval_cut & H3K4me3$V4>0)] = paste0("UP-",length(which(H3K4me3$V5<pval_cut & H3K4me3$V4>0)))

  colnames(all)[8] = "ann"
  
  # mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(16)
  
  g = ggplot(all) + geom_col(aes(x=sig,y=1,fill=ann),position="fill") +
    theme_bw()  +labs(title=i)+scale_fill_manual(values =brewer.pal(7, "Set2"))
  #print(g)
  p[[i]] = g
}

setwd("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/overlap_wi_histone")
pdf(paste0(tissue,".custom_ann_all_together.1p.pdf"), height = 18, width = 20)
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]])
dev.off()
