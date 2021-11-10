#setwd("/mnt/silencer2/home/lamaral/lamaral/projects/Aging/histone_overlap/")

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

  dd$Clade = factor(dd$Clade, levels = c("N/A", "Other" ,"Myeloid","Lymphoid", "Muscle", "ExN" , "InN", "Glia"),ordered = T)

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
#  sig_label[which(sig<=0.025)] = "*"
  sig_label[which(sig<=0.001)] = "**"
  sig_label[which(sig<=0.00001)] = "***"
  sig_label[which(sig>0.025)] = ""
#  p <- ggplot(dd, aes(x = OR, y = id, color = Clade)) +
  p <- ggplot(dd, aes(x = -log10(p), y = id, color = Clade,size=OR)) +
    geom_vline(aes(xintercept = 2), size = .25, linetype = "dashed") + scale_color_brewer(palette="Dark2")+
#    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .2, color = "gray50") +
    geom_point() +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    #    scale_y_continuous(breaks = seq(1,nrow(dd),1), labels = dd$id) +
    scale_x_continuous() + coord_trans(x = "sqrt") +
#    ylab("") + xlim(c(xmin , max))+
    xlab("-log10 P-value") +
    annotate(geom = "text", y =.6, x = xmax-1, label ="", size = 3.5, hjust = 0) +
    ggtitle(paste(mark, dir, "Overlap Enrichment")) + theme(text = element_text(face="bold",  size=12, angle=0),axis.text.y = element_text(face="bold",  size=12, angle=0),
                                                            axis.text.x = element_text(face="bold",  size=12, angle=0)) +
    facet_grid(Clade~.,scales="free_y",space="free_y")
#    geom_text(aes(label=sig_label), vjust=-0.1, hjust = -0.1) +  guides(colour = guide_legend(reverse=T))
  p
  return(p)
}

tab = read.csv(paste("OR_tables/", "H3K9me3","_OR_tab.txt" , sep = ""), header = T, sep = "\t")
meta = read.delim("../../celltype_annotation.txt")
tab$Clade = meta$Clade[match(paste(tab$Tissue,tab$Name),paste(meta$Tissue,meta$Name))]
#tab = tab[which(tab$up_overlap>=5),]

up = getORBoxPlotLog(tab,"Up","H3K9me3")
pdf("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/overlap_wi_histone/all_tissues/test.H3K9me3.up.all_tissues.pdf", height = 10, width = 8)
print(up)
dev.off()




