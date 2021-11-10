
## requires loj_file : file containing all peaks in a cell type, p-values and logFC across aging, and
## whether it overlaps with a mark (generated with bedtools intersect -loj)

# loj files can be found at /mnt/tscc/lamaral/projects/Aging/[tissue]/overlap/peaks.bed/

### get the overlap enrichment for a single cluster
### computes fischer's exact p-value and OR for up vs all peaks, down vs all peaks enrichment
getSig <- function(loj_file,cl, tissue, p_val) {
  cat(cl,"\n")
  loj = read.csv(loj_file, sep = "\t", header = F)
  loj$id = paste0(loj$V1,":",loj$V2,"-",loj$V3, sep = "")
  ## remove duplicated ATAC peaks 
  if (length(which(duplicated(loj$id)))>0) {
    loj = loj[-which(duplicated(loj$id)), ]
  }
  
### to remove dup ATAC peaks in same histone peak
  #table(loj$V10)[order(table(loj$V10), decreasing = T)[1:10]]
  #  if(length(which(duplicated(loj$V10) & loj$V10 != "."))>0) {
  #    loj = loj[-which(duplicated(loj$V10) & loj$V10 != "."),]
  #  }
  #cat(dim(loj), "\n")
  
  loj = loj[-grep("chrM", loj$V1),]
  diff = loj
  tfbound = as.vector(loj$V7)
  tfbound[tfbound != "."] <- paste("Bound Peaks")
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
  
  fisher.test(dwn_all)
  fdn = fisher.test(dwn_all)
  
  cat(fup$estimate,fdn$estimate, "\n")
  cat(fup$p.value,fdn$p.value, "\n-----\n")
  
  up_info = c(cl,"Up",fup$p.value,fup$estimate,fup$conf.int)
  down_info = c(cl, "Down",fdn$p.value,fdn$estimate,fdn$conf.int)
  ret_df = data.frame(rbind(up_info, down_info))
  ret_df$up_total = up_Age_total
  ret_df$down_total = dwn_Age_total
  ret_df$up_overlap = up_Age_tf_Bound
  ret_df$down_overlap = dwn_Age_tf_Bound
  colnames(ret_df) = c("cl","dir","p","OR", "CILow", "CIHigh", "up_total", "down_total", "up_overlap", "down_overlap")
  return(ret_df)
}

# get the overlap enrichment for all clusters by calling getSig
getSigTable <- function(tissue, mark,p_val) {
  setwd(paste("/mnt/tscc/lamaral/projects/Aging/",tissue,"/overlap/", sep = ""))
  final_table=data.frame()
  final_table  = rbind(final_table,c("cl","dir","p","OR", "CILow", "CIHigh"))
  colnames(final_table) = c("cl","dir","p","OR", "CILow", "CIHigh")
  final_table = final_table[-1,]
  files = list.files(path="peaks.bed",pattern=paste(mark,".loj.bed", sep = ""),full.names=T)
  for (cl in 1:length(files)) {
    file = paste("peaks.bed/",cl,".",mark,".loj.bed", sep = "")
    curr = getSig(file,cl, tissue,p_val)
    final_table = rbind(final_table , curr)
  }
  if (tissue == "FC") {
    group = c("Neuron", "Neuron", "Glial","Glial", "Neuron","Neuron","Neuron",
              "Neuron","Neuron","Neuron","Glial","Neuron","Neuron","Glial",
              "Neuron","Glial", "Glial", "Neuron")
    id = c("Ex.L2-3.Cux2", "Ex.L2-3.Cux1", "Olg.Mog", "Ast.Apoe", "Inh.Npy", 
           "Ex.L5.Capsl", "Ex.L6.Sulf1.Syt6", "Inh.Meis2", "Inh.Sst", 
           "Inh.Six3","Micg.C1qb", "Ex.L6.Sulf1.Nnat", "Ex.L5.Kcnn2",
           "OPC.Pdgfra", "Ex.L5.Tszh2", "Peri.Foxd1", "Endo.Tie1", "Inh.Lmc38")
  } else if (tissue == "DH") {
    group = c("Neuron", "Glial", "Neuron", "Neuron", "Neuron", "Neuron", 
              "Neuron", "Neuron", "Glial","Glial","Glial","Glial","Glial","Glial")
    id = c("DG.Prox1.Penk", "Olg.Mog", "CA1.Fibdc1", "DG.Prox1.Card6", "DG.Prox1.Card6.2", "Sub_Ent", "Inh.Gad2",
           "CA2/3.Cldn22", "Ast.Apoe", "Mcg.C1qb", "Chor.Plex.Ttr","OPC.Pdgfra", "Endo.Tie1","Peri.Foxd1")
  } else if (tissue == "LM") {
    group = c("Skeletal muscle", "Skeletal muscle", "Skeletal muscle", "Skeletal muscle", "Fibroblast", 
              "Skeletal muscle", "Kerat","Endothelial","Imm/blood","Smooth muscle","Fat", "Muscle stem cell", "Imm/blood", 
              "Skeletal muscle", "Skeletal muscle", "Imm/blood", "Imm/blood")
    id = c("Skm.Fhl3.Kcnf1", "Skm.Fhl3.Myh1.Myl3", "Skm.Fhl3.LowQuality", "Skm.Fhl3.Kcnf1.1", "FB.Dcn", 
           "Skm.Fhl3.Atp1b4", "Kerat.Krt5","Endo.Tie1","Macg.C1qb","SMC.Myh11","Fat.Cidec", "MSC.Pax7",
           "Bcell.Cd79a", "Skm.Fhl3.Crym.Ache", "Skm.Fhl3.Wdr31", "Neutrophils.S100a8", "Eryth.Hbb-bs")
  }
  else if (tissue == "HT") {
    group = c("Cardiomyocyte", "Fibroblast", "Endothelial", "Endothelial", "Imm/blood", "Smooth muscle", "Endothelial","Imm/blood","Fat","Imm/blood","Cardiomyocyte")
    id = c("CM.Myh6", "FB.Sox9.Dcn", "Endo.Tie1", "Low_quality?", "Macg.Mgl2", "SM.Myh11", "Endo.Tie1.1",
           "Bcell.Cd79a","Fat.Cidec","Tcell.Cd96","CM.FB.doublet?")
  }
  final_table$group = rep(group,each = 2)
  final_table$id = rep(id,each = 2)
  return(final_table)
}

# plots one odds ratio plot
getORBoxPlotLog <- function(OR_table, dir, tissue, mark, xmin,xmax) {
  dd = data.frame(OR_table[which(OR_table$dir==dir),],stringsAsFactors = F)
  dd$cl=as.numeric(levels(dd$cl))[dd$cl]
  dd$p=as.numeric(levels(dd$p))[dd$p]
  dd$p = dd$p * nrow(dd)
  dd$OR=(as.numeric(levels(dd$OR))[dd$OR])
  dd$CILow=(as.numeric(levels(dd$CILow))[dd$CILow])
  dd$CIHigh=(as.numeric(levels(dd$CIHigh))[dd$CIHigh])
  dd$diff = log10(dd$CIHigh+.01)-log10(dd$CILow+.01)
  if(length(which(dd$diff > 2 & (sign(log(dd$CILow+.01) * log(dd$CIHigh)) == -1)))>0) {
    dd = dd[-which(dd$diff > 2 & (sign(log(dd$CILow+.01) * log(dd$CIHigh)) == -1)),]    
  }
  if(length(which(dd$CIHigh>xmax))>0){
    dd[which(dd$CIHigh>xmax),"CIHigh"] = xmax
  }
  if(length(which(dd$CILow<xmin))>0){
    dd[which(dd$CILow<xmin),"CILow"] = xmin
  }
  if(length(which(dd$OR<xmin))>0){
    dd[which(dd$OR<xmin),"OR"] = xmin
  }
  if(length(which(dd$OR>xmax))>0){
    dd[which(dd$OR>xmax),"OR"] = xmax
  }
  sig = signif(dd$p,2)
  sig_label = sig
  sig_label[which(sig<=0.05)] = "*"
  sig_label[which(sig<=0.001)] = "**"
  sig_label[which(sig<=0.00001)] = "***"
  sig_label[which(sig>0.05)] = ""
  if (tissue == "FC" || tissue == "DH") {
    dd$id = factor(dd$id, levels = dd$id[order(dd$group)])
  }
  else if (dir == "Up") {
    dd$id = factor(dd$id, levels = dd$id[order(dd$group, dd$up_total+dd$down_total, decreasing = T)])
  } else {
    dd$id = factor(dd$id, levels = dd$id[order(dd$group, dd$up_total+dd$down_total, decreasing = T)])
  }
  p <- ggplot(dd, aes(x = OR, y = id, color = group)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = CIHigh, xmin = CILow), size = .5, height = .2, color = "gray50") +
    geom_point(size = 3.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    #    scale_y_continuous(breaks = seq(1,nrow(dd),1), labels = dd$id) +
    scale_x_continuous(breaks = 10) + coord_trans(x = "log10") +
    ylab("") + xlim(c(xmin , xmax))+
    xlab("Odds ratio") +
    annotate(geom = "text", y =.6, x = xmax-1, label ="", size = 3.5, hjust = 0) + 
    ggtitle(paste(tissue, mark, dir)) + theme(text = element_text(face="bold",  size=12, angle=0),axis.text.y = element_text(face="bold",  size=12, angle=0),
                                              axis.text.x = element_text(face="bold",  size=12, angle=0))+ 
    geom_text(aes(label=sig_label), vjust=-0.1, hjust = -0.1)
  p
  return(p)
}

# plots 4 panel odds ratio plot (FC,DH up and down)
# xmax is the maximum odds ratio on the x axis
get4PanelRes <- function(mark,xmax,p_val) {
  FC_table = getSigTable("FC", mark,p_val)
  DH_table = getSigTable("DH", mark,p_val)
  if (length(which(DH_table$up_total < 10 | DH_table$down_total < 10)>1)){
    DH_table = DH_table[-which(DH_table$up_total < 10 | DH_table$down_total < 10), ]
  }
  if (length(which(FC_table$up_total < 10 | FC_table$down_total < 10)>1)){
    FC_table = FC_table[-which(FC_table$up_total < 10 | FC_table$down_total < 10), ]
  }
  fd = getORBoxPlotLog(FC_table, "Down", "FC", mark,10^-(log10(xmax)),xmax)
  dd = getORBoxPlotLog(DH_table, "Down", "DH", mark,10^-(log10(xmax)),xmax)
  fu = getORBoxPlotLog(FC_table, "Up", "FC", mark,10^-(log10(xmax)),xmax)
  du = getORBoxPlotLog(DH_table, "Up", "DH", mark,10^-(log10(xmax)),xmax)
  print(grid.arrange(fd,dd, fu,du))
  #return(c(fd,dd, fu,du))
  return(list(FC_table, DH_table))
}


# code just for generating the H3K9me3 reprocessed plots in FC, DH
setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/histone_mark_overlap/H3K9me3_reprocessed/")
pdf("Odds_ratio_plot_H3K9me3_before_and_after_reprocessing_correct.pdf", height = 10, width = 12)
H3K9me3 = get4PanelRes("H3K9me3", 25, 0.01)
H3K9me3_reprocessed = get4PanelRes("H3K9me3_reprocessed", 25, 0.01)
dev.off()
write.csv(H3K9me3[[1]], "FC_H3K9me3_overlap_table.txt",quote = F, row.names = F)
write.csv(H3K9me3[[2]], "DH_H3K9me3_overlap_table.txt",quote = F, row.names = F)
write.csv(H3K9me3_reprocessed[[1]], "FC_H3K9me3_reprocessed_overlap_table.txt", quote = F, row.names = F)
write.csv(H3K9me3_reprocessed[[2]], "DH_H3K9me3_reprocessed_overlap_table.txt", quote = F, row.names = F)



# original code for calling all of the marks, 

pdf("/mnt/tscc/lamaral/projects/Aging/aging_share/overlap/Odds_ratio_plots_FC_DH.pdf", height = 10, width = 12)
plot_CTCF = get4PanelRes("CTCF",6,0.01)
plot_H3K9me3 = get4PanelRes("H3K9me3",14,0.01)
plot_H3K27me3 = get4PanelRes("H3K27me3",6,0.01)
plot_H3K4me3 = get4PanelRes("H3K4me3",8,0.01)
plot_H3K4me1 = get4PanelRes("H3K4me1",6,0.01)
plot_H3K36me3 = get4PanelRes("H3K36me3",6,0.01)
plot_H3K27ac = get4PanelRes("H3K27ac",8,0.01)
plot_H3K9ac = get4PanelRes("H3K9ac",8,0.01)
plot_satellite = get4PanelRes("mandy_blacklist", .01, 14,0.01)
dev.off()



# for plotting LM and HT 4 panel plot
get4PanelResLMHT <- function(mark,xmax, p_val) {
  LM_table = getSigTable("LM", mark,p_val)
  HT_table = getSigTable("HT", mark,p_val)
  LM_table = LM_table[-which(LM_table$up_total < 10 | LM_table$down_total < 10), ]
  HT_table = HT_table[-which(HT_table$up_total < 10 | HT_table$down_total < 10), ]
  fd = getORBoxPlotLog(LM_table, "Down", "LM", mark,10^-(log10(xmax)),xmax)
  dd = getORBoxPlotLog(HT_table, "Down", "HT", mark,10^-(log10(xmax)),xmax)
  fu = getORBoxPlotLog(LM_table, "Up", "LM", mark,10^-(log10(xmax)),xmax)
  du = getORBoxPlotLog(HT_table, "Up", "HT", mark,10^-(log10(xmax)),xmax)
  print(grid.arrange(fd,dd, fu,du))
  return(c(fd,dd, fu,du))
}


pdf("/mnt/tscc/lamaral/projects/Aging/aging_share/overlap/Odds_ratio_plots_HT_LM.pdf", height = 10, width = 12)
plot_CTCF = get4PanelResLMHT("CTCF",6,0.01)
plot_H3K9me3 = get4PanelResLMHT("H3K9me3",15,0.01)
plot_H3K27me3 = get4PanelResLMHT("H3K27me3",8,0.01)
plot_H3K4me3 = get4PanelResLMHT("H3K4me3",6,0.01)
plot_H3K4me1 = get4PanelResLMHT("H3K4me1",6,0.01)
plot_H3K36me3 = get4PanelResLMHT("H3K36me3",6,0.01)
plot_H3K27ac = get4PanelResLMHT("H3K27ac",8,0.01)
plot_H3K9ac = get4PanelResLMHT("H3K9ac",6,0.01)
dev.off()

