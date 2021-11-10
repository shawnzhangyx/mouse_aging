################## count the number of differential peaks using the same cut-off.
library(pheatmap)
library(reshape2)
library(ggplot2)
library(GenomicRanges)

# gets differential peak data for all cell types in tissue in one table
getDiffTable = function(tissue) {
  setwd(paste("/mnt/tscc/yanxiao/projects/mouse_aging/analysis/snapATAC/",tissue,"/age_diff_edgeR.snap",sep = ""))
  b = read.delim("../../../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
  id = b$Name[which(b$Tissue==tissue)]

  ## read the meta info. 
  meta = read.delim(paste("../",tissue,".pool.barcode.meta_info.txt",sep = ""))
  
  files = list.files(pattern=".edger.txt")
  max_cluster = length(files)
  dat = list()
  
  
  for (cl in 1:max_cluster) {
    dat[[cl]] = data.frame(fread(paste0(cl,".edger.txt")))
    dat[[cl]]$cluster = cl
    # Select the top 1% of peaks as diff. 
    dat[[cl]] = dat[[cl]][which(dat[[cl]]$PValue < quantile(dat[[cl]]$PValue,0.01)),]
  }
  
  dat2 = do.call(rbind, dat)
  # Number of Shared differential peaks.
   diff = dat2  
  # diff = dat2[dat2$PValue < 0.01,]
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
HT_diff = getDiffTable("HT")
colnames(HT_diff) = gsub("HT_", "", colnames(HT_diff))
LM_diff = getDiffTable("LM")
colnames(LM_diff) = gsub("LM_", "", colnames(LM_diff))
BM_diff = getDiffTable("BM")
colnames(BM_diff) = gsub("BM_", "", colnames(BM_diff))


#all tissues diff table
diff = rbind(DH_diff,FC_diff, LM_diff, HT_diff,BM_diff)

# shared diff bar plot
#shared_diff = data.frame(table(table(diff$X))[-1])
#setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/diff_peak_overlap/all_tissues_together/")
#pdf("All_tissues_Diff_peaks_shared_by_clusters.001.pdf")
#ggplot(shared_diff) + geom_col(aes(x=Var1,y=Freq)) +
#  geom_text(aes(x=Var1,y=Freq,label=Freq),nudge_y=100)+
#  xlab("# of Cell Types") + ylab("# of diff peaks")
#dev.off()

diff$tissue_celltype = paste0(diff$tissue,"_",diff$cell_type, sep="")
diff$tissue_cluster = paste0(diff$tissue,"_",diff$cluster, sep="")
meta = read.delim("../../../../aging_share/figures/celltype_annotation.txt",stringsAsFactors=F)
diff = diff[-which(diff$tissue_celltype %in% paste0(meta$Tissue,"_",meta$Name)[is.na(meta$Clade)]),]



peak_cnt = data.frame(table(diff$tissue_celltype))
peak_cnt = peak_cnt[order(-peak_cnt$Freq),]

t2 = table(diff$PValue < 1, diff$logFC>0, diff$tissue_celltype, diff$tissue_cluster, diff$tissue, diff$cluster)
mt2 =melt(t2)

# # diff peak per cell type/tissue bar plot
setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")

pdf("all_tissues_number_of_diff_peaks_1pc.pdf",height = 12,width = 10)
ggplot(mt2) +
#  geom_col(aes(x=factor(Var3, levels = rev.Vector(unique(diff$tissue_celltype))),
geom_col(aes(x=factor(Var3, levels = rev(peak_cnt$Var1)),
      y=value,fill=Var2),position="dodge") +
  xlab("Cluster") +
  scale_fill_discrete(name="Up in Aging") +
  coord_flip()
dev.off()


mat = matrix(0,nrow=length(unique(diff$tissue_celltype)),ncol=length(unique(diff$tissue_celltype)))
x = 1
for (i in unique(diff$tissue_celltype)){
  y=1
  set1 = diff[diff$tissue_celltype==i,"V1"] 
  gr1 = GRanges(set1)
  for (j in unique(diff$tissue_celltype)){
    print(paste(i,j))
    set2 = diff[diff$tissue_celltype==j,"V1"]
    gr2 = GRanges(set2) 
    mat[x,y] = length(intersect(gr1,gr2))/length(union(gr1,gr2))
    #mat[x,y] = length(intersect(set1,set2))/length(union(set1,set2))
    y = y+1
  }
  x = x+1
}

#melted = melt(mat)
colnames(mat) = unique(diff$tissue_celltype)
rownames(mat) = colnames(mat)

write.csv(mat,"All_tissues.diffpeak_share.table.1pc.csv")
#mat = read.csv("All_tissues.diffpeak_share.table.1pc.csv",row.names=1)
mat[which(mat==1, arr.ind = T)] = 0.02 
mat[which(mat>0.02, arr.ind = T)] = 0.02
#logmat = log10(mat+0.0001)

anno_row = data.frame(tissue=substr(colnames(mat),1,2))
rownames(anno_row) = colnames(mat)
anno_col = data.frame(clade=meta$Clade[match(colnames(mat),paste0(meta$Tissue,"_",meta$Name))])
rownames(anno_col) = colnames(mat)

# jaccard similarity for diff peaks bar plot
setwd("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/diffpeak_num_overlap/")
pdf("all_tissues_jaccard_heat_1pc.pdf", height = 14, width = 16)
#pheatmap(logmat, main = paste( "All tissues Jaccard of differential peaks \n log10(jaccard+0.0001)"),legend_labels = "log10(jaccard +0.0001)")
#pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks \n max 0.02"),annotation_row =anno_row,annotation_col = anno_col)
pheatmap(mat, main = paste( "All tissues Jaccard of differential peaks \n max 0.02"),
  clustering_distance_rows = "correlation", 
  clustering_distance_cols = "correlation",
  border_color="black",
  annotation_row =anno_row,
  annotation_col = anno_col,
  )
dev.off()




