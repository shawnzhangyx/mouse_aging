################## count the number of differential peaks using the same cut-off.
library(pheatmap)
library(reshape2)
library(ggplot2)
library(GenomicRanges)

# gets differential peak data for all cell types in tissue in one table
getDiffTable = function(tissue) {
  setwd(paste("/mnt/tscc/yanxiao/projects/mouse_aging/analysis/snapATAC/",tissue,"/age_diff_edgeR.snap",sep = ""))
  if (tissue == "DH") {
    id = c("DG.Prox1.Penk", "Olg.Mog", "CA1.Fibdc1", "DG.Prox1.Card6", "DG.Prox1.Card6.2", "Sub_Ent", "Inh.Gad2",
           "CA2/3.Cldn22", "Ast.Apoe", "Mcg.C1qb", "Chor.Plex.Ttr","OPC.Pdgfra", "Endo.Tie1","Peri.Foxd1")
  }
  if (tissue == "FC") {
    id = c("Ex.L2-3.Cux2", "Ex.L4.Cux1", "Olg.Mog", "Ast.Apoe", "Inh.Npy", 
           "Ex.L5.Capsl", "Ex.L6.Sulf1.Syt6", "Inh.Meis2", "Inh.Sst", 
           "Inh.Six3","Micg.C1qb", "Ex.L6.Sulf1.Nnat", "Ex.L5.Kcnn2",
           "OPC.Pdgfra", "Ex.L5.Tszh2", "Peri.Foxd1", "Endo.Tie1", "Inh.Lmc38")
  }
  
  ## read the meta info. 
  meta = read.delim(paste("../",tissue,".pool.barcode.meta_info.txt",sep = ""))
  
  files = list.files(pattern=".edger.txt")
  max_cluster = length(files)
  dat = list()
  
  
  for (cl in 1:max_cluster) {
    dat[[cl]] = read.csv(paste0(cl,".edger.txt"))
    dat[[cl]]$cluster = cl
  }
  
  dat2 = do.call(rbind, dat)
  # Number of Shared differential peaks.
  diff = dat2[dat2$PValue < 0.001,]
  diff$tissue = tissue
  diff$cell_type = ""
  for (i in 1:nrow(diff)) {
    diff$cell_type[i] = id[diff$cluster[i]]
  }
  return(diff)
}
DH_diff = getDiffTable("DH") 
colnames(DH_diff) = gsub("DH_", "", colnames(DH_diff))
FC_diff = getDiffTable("FC")
colnames(FC_diff) = gsub("FC_", "", colnames(FC_diff))

#all tissues diff table
diff = rbind(DH_diff,FC_diff)


diff$tissue_celltype = paste0(diff$tissue,"_",diff$cell_type, sep="")
diff$tissue_cluster = paste0(diff$tissue,"_",diff$cluster, sep="")


mat = matrix(0,nrow=length(unique(diff$tissue_celltype)),ncol=length(unique(diff$tissue_celltype)))
x = 1
for (i in unique(diff$tissue_celltype)){
  y=1
  set1 = diff[diff$tissue_celltype==i,"X"] 
  gr1 = GRanges(set1)
  for (j in unique(diff$tissue_celltype)){
    print(paste(i,j))
    set2 = diff[diff$tissue_celltype==j,"X"]
    gr2 = GRanges(set2) 
    mat[x,y] = length(intersect(gr1,gr2))/length(union(gr1,gr2))
    #mat[x,y] = length(intersect(set1,set2))/length(union(set1,set2))
    y = y+1
  }
  x = x+1
}

melted = melt(mat)
mat2 = mat

mat2[which(mat2==1, arr.ind = T)] = NA
mat2[which(mat2>0.05, arr.ind = T)] = 0.05
rownames(mat2) = unique(diff$tissue_celltype)
colnames(mat2) = unique(diff$tissue_celltype)
logmat2 = log10(mat2+0.00001)

# jaccard similarity for diff peaks bar plot
setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")
pdf("DH_FC_jaccard_heat_ph.001.pdf", height = 12, width = 12)
pheatmap(logmat2, main = paste( "All tissues Jaccard of differential peaks \n log10(jaccard+0.0001)"),legend_labels = "log10(jaccard +0.0001)")
pheatmap(mat2, main = paste( "All tissues Jaccard of differential peaks \n max 0.05"))
dev.off()

