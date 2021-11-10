# to create odds ratio plots from OR tables. 
setwd("/mnt/silencer2/home/lamaral/lamaral/projects/Aging/histone_overlap/")

getORBoxPlotLog <- function(OR_table, dir, mark) {
  dd = data.frame(OR_table[which(OR_table$dir==dir),],stringsAsFactors = F)
  # I uncomment this to make the brain only plots 
  #dd = dd[which(dd$Tissue %in% c("DH","FC")),]
  #dd$Clade[which(dd$Clade=="Immune")] = "Glia"
  
  dd$id = paste0(dd$Tissue,".",dd$Name)
  dd$id = make.unique(dd$id)
  dd$Clade = as.character(dd$Clade)
  dd$Clade[which(is.na(dd$Clade))] = "N/A"
  dd$Clade[which(dd$Clade=="")] = "Other"
  
  dd$Clade = factor(dd$Clade, levels = c("N/A", "Other" ,  "Immune", "Muscle", "ExN" , "InN", "Glia"),ordered = T)
  
  dd = dd[order(dd$Clade, dd$OR,na.last = F),]
  dd$id  = factor(dd$id, levels = dd$id)
  
  xmax = sort(dd$CIHigh,decreasing = T)[which(sort(dd$CIHigh,decreasing = T) < 100)[1]]
  xmin = sort(dd$CILow)[which(sort(dd$CILow) >0.01)[1]]
  
  xmax = max(xmax, 1/xmin)
  xmin = 1/xmax
  
  if (mark == 'H3K9me3_reprocessed') {
    mark = "H3K9me3"
  }
  if(length(which(dd$CILow<xmin | dd$CIHigh>xmax))>0) {
    dd =dd[-which(dd$CILow<xmin | dd$CIHigh>xmax),]
  }
  dd = dd[-which(dd$id %in% c("HT.doublet","HT.Low_quality","LM.Low_quality","FC.Inh.Lrrc38.Plcxd3.Prok2","LM.Eryth")),]
  
  sig = dd$p
  sig_label = sig
  sig_label[which(sig<=0.025)] = "*"
  sig_label[which(sig<=0.001)] = "**"
  sig_label[which(sig<=0.00001)] = "***"
  sig_label[which(sig>0.025)] = ""
  
  p <- ggplot(dd, aes(x = OR, y = id, color = Clade)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + scale_color_brewer(palette="Dark2")+
    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .2, color = "gray50") +
    geom_point(size = 3.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    #    scale_y_continuous(breaks = seq(1,nrow(dd),1), labels = dd$id) +
    scale_x_continuous(breaks = 10) + coord_trans(x = "log10") +
    ylab("") + xlim(c(xmin , xmax))+
    xlab("Odds ratio") +
    annotate(geom = "text", y =.6, x = xmax-1, label ="", size = 3.5, hjust = 0) + 
    ggtitle(paste(mark, dir, "Overlap Enrichment")) + theme(text = element_text(face="bold",  size=12, angle=0),axis.text.y = element_text(face="bold",  size=12, angle=0),
                                                            axis.text.x = element_text(face="bold",  size=12, angle=0))+ 
    geom_text(aes(label=sig_label), vjust=-0.1, hjust = -0.1) +  guides(colour = guide_legend(reverse=T))
  p
  return(p)
}

marks = c("H3K9me3","CTCF","H3K4me3", "H3K27ac","H3K4me1", "H3K4me2", "H3K36me3","H3K9ac","H3K27me3",  "A_comp", "B_comp")
for(mark in marks) {
  tab = read.csv(paste("OR_tables/", mark,"_OR_tab.txt" , sep = ""), header = T, sep = "\t")
  up = getORBoxPlotLog(tab,"Up",mark)
  down = getORBoxPlotLog(tab,"Down",mark)
  #pdf(paste("brain_plots/",mark,"_all_tissues.pdf",sep =""), height = 12, width = 8)
  pdf(paste("plots/",mark,"_all_tissues.pdf",sep =""), height = 14, width = 8)
  print(up)
  print(down)
  dev.off()
}
  
