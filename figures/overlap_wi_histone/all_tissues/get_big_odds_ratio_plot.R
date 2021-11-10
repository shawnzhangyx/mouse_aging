

setwd("/mnt/tscc/lamaral/projects/Aging/histone_overlap/")
## requires overlap file : file containing all peaks in a cell type, p-values and logFC across aging, and
## whether it overlaps with a mark (generated with bedtools intersect -c)

# over files can be found at /mnt/tscc/lamaral/projects/Aging/histone_marks/[tissue]/[tissue]_pverlap

key = read.csv("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/celltype_annotation.txt", sep = "")


### get the overlap enrichment for a single cluster
### computes fischer's exact p-value and OR for up vs all peaks, down vs all peaks enrichment
getSig <- function(file) {
  over = read.csv(file, sep = "\t", header = F)
  chr = over$V1
  start = over$V2
  end = over$V3
  over$id = paste0(over$V1,":",over$V2,"-",over$V3, sep = "")
  
  over = over[-grep("chrM", over$V1),]
  diff = over
  tfbound = as.vector(over$V7)
  tfbound[tfbound != "0"] <- paste("Bound Peaks")
  tfbound[tfbound != paste("Bound Peaks")] <- "ALL PEAKS"
  diff$tf_bound = as.factor(tfbound)
  colnames(diff)[4] <- "Log2.Fold.Change"
  colnames(diff)[5] <- "p.value"
  diff = diff[order(diff$p.value),]
  up_Age_tf_Bound = nrow(diff[which(
    diff$tf_bound == paste("Bound Peaks") &
      diff$p.value < diff[nrow(diff)*0.01,"p.value"] &
      diff$Log2.Fold.Change > 0
  ), ])
  
  up_Age_total =  nrow(diff[which(
    diff$p.value < diff[nrow(diff)*0.01,"p.value"] &
      diff$Log2.Fold.Change > 0
  ), ])
  
  dwn_Age_tf_Bound = nrow(diff[which(
    diff$tf_bound == paste("Bound Peaks") &
      diff$p.value < diff[nrow(diff)*0.01,"p.value"] &
      diff$Log2.Fold.Change < 0
  ), ])
  
  dwn_Age_total =  nrow(diff[which(
    diff$p.value < diff[nrow(diff)*0.01,"p.value"] &
      diff$Log2.Fold.Change < 0
  ), ])
  
  all_tf_Bound = nrow(diff[which(diff$tf_bound == paste("Bound Peaks")), ])
  
  total_peaks = nrow(diff)
  
  up_all = rbind(
    c(up_Age_tf_Bound, up_Age_total - up_Age_tf_Bound),
    c(all_tf_Bound, total_peaks - all_tf_Bound)
  )
  
  fisher.test(up_all,conf.int = T)
  fup = fisher.test(up_all)
  
  dwn_all = rbind(
    c(dwn_Age_tf_Bound, dwn_Age_total - dwn_Age_tf_Bound),
    c(all_tf_Bound, total_peaks - all_tf_Bound)
  )
  
  fisher.test(dwn_all,conf.int = T)
  fdn = fisher.test(dwn_all)
  
  cat(fup$estimate,fdn$estimate, "\n")
  cat(fup$p.value,fdn$p.value, "\n-----\n")
  
  up_info = c("Up",fup$p.value,fup$estimate,fup$conf.int)
  down_info = c("Down",fdn$p.value,fdn$estimate,fdn$conf.int)
  ret_df = data.frame(rbind(up_info, down_info))
  ret_df$up_total = up_Age_total
  ret_df$down_total = dwn_Age_total
  ret_df$up_overlap = up_Age_tf_Bound
  ret_df$down_overlap = dwn_Age_tf_Bound
  ret_df$all_total = total_peaks
  ret_df$all_overlap = all_tf_Bound
  colnames(ret_df) = c("dir","p","OR", "CILow", "CIHigh", "up_total", "down_total", "up_overlap", "down_overlap","all_total","overlap_total")
  return(ret_df)
}

# get the overlap enrichment for all clusters by calling getSig
getSigTable <- function(mark) {
  final_table=data.frame()
  final_table  = rbind(final_table,c("dir","p","OR", "CILow", "CIHigh","Tissue", "Cluster", "Name","Clade","file"))
  colnames(final_table) = c("dir","p","OR", "CILow", "CIHigh")
  final_table = final_table[-1,]
  for (i in 1:nrow(key)) {
    setwd(paste("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/overlap_wi_histone/all_tissues/",key$Tissue[i],"_overlap/", sep = ""))
    cat(paste("intersect_",key$Cluster[i],"_",mark, sep = ""),"\n")
    file = list.files(path = ".", pattern = paste("intersect_",key$Cluster[i],"_",mark, sep = ""))[1]
    if(is.na(file)) {
      next()
    }
    curr = getSig(file)
    curr = cbind(curr,key$Tissue[i], key$Cluster[i], key$Name[i],key$Clade[i])
    colnames(curr)[12:15] = c("Tissue", "Cluster", "Name","Clade")
    curr$file = file
    final_table = rbind(final_table , curr)
  }
  return(final_table)
}

# plots one odds ratio plot
getORBoxPlotLog <- function(OR_table, dir, mark) {
  dd = data.frame(OR_table[which(OR_table$dir==dir),],stringsAsFactors = F)
  # I uncomment this to make the brain only plots 
  #dd = dd[which(dd$Tissue %in% c("DH","FC")),]
  #dd$Clade[which(dd$Clade=="Immune")] = "Glia"
  dd$p=as.numeric(levels(dd$p))[dd$p]
  dd$OR=(as.numeric(levels(dd$OR))[dd$OR])
  dd$CILow=(as.numeric(levels(dd$CILow))[dd$CILow])
  dd$CIHigh=(as.numeric(levels(dd$CIHigh))[dd$CIHigh])
  dd$id = paste0(dd$Tissue,".",dd$Name)
  dd$id = make.unique(dd$id)
  dd$Clade = as.character(dd$Clade)
  dd$Clade[which(is.na(dd$Clade))] = "N/A"
  dd$Clade[which(dd$Clade=="")] = "Other"
  
  dd$Clade = factor(dd$Clade, levels = c("N/A", "Other" ,  "Immune", "Muscle", "ExN" , "InN", "Glia"),ordered = T)
  
  dd = dd[order(dd$Clade, dd$OR,na.last = F),]
  dd$id  = factor(dd$id, levels = dd$id)
  sig = dd$p
  sig_label = sig
  sig_label[which(sig<=0.025)] = "*"
  sig_label[which(sig<=0.001)] = "**"
  sig_label[which(sig<=0.00001)] = "***"
  sig_label[which(sig>0.025)] = ""
  
  xmax = sort(dd$CIHigh,decreasing = T)[which(sort(dd$CIHigh,decreasing = T) != Inf)[1]]
  xmin = sort(dd$CILow)[which(sort(dd$CILow) != 0)[1]]

  xmax = max(xmax, 1/xmin)
  xmin = 1/xmax
  
  p <- ggplot(dd, aes(x = OR, y = id, color = Clade)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + scale_color_brewer(palette="Dark2")+
    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .2, color = "gray50") +
    geom_point(size = 3.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
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
  tab = getSigTable(mark)
  setwd("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/overlap_wi_histone/all_tissues/OR_tables/")
  write.table(tab, file = paste(mark, "_OR_tab.txt", sep = ""), sep = "\t", quote = F, row.names = F)
  setwd("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/overlap_wi_histone/all_tissues/plots/")
  up = getORBoxPlotLog(tab,"Up",mark)
  dp = getORBoxPlotLog(tab,"Down",mark)
  pdf(paste(mark,"_all_tissues.pdf",sep =""), height = 10, width = 8)
  print(up)
  print(dp)
  dev.off()
}
