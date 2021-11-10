library(GenomicRanges)
library(pheatmap)
library(matrixStats)
library(RColorBrewer)

key = read.csv("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/celltype_annotation.txt", sep = "",stringsAsFactors=F)
key = key[!is.na(key$Clade),]

setwd("/mnt/tscc/lamaral/projects/Aging/GWAS/")
# lift is all of the liftover SNPs
#lift = read.csv("bed_LD_snp_files_wtag/hglft_genome_wtag.bed", sep = "\t", header = F,stringsAsFactors = F)
lift = read.csv("old/hglft_genome_all_LD_snps.bed", sep = "\t", header = F,stringsAsFactors = F)

traits_to_remove = names(which(table(lift$V4)<30))
traits_to_keep = read.delim("/projects/ps-renlab/yanxiao/projects/mouse_aging/data/all_aging_traits_EFO.curated.txt",header=F,stringsAsFactors=F)$V1
traits_to_keep = gsub(" ",".",traits_to_keep)
traits_to_keep = gsub("-",".",traits_to_keep)
traits_to_keep = gsub("'",".",traits_to_keep)


countOccurance <- function(text, vect) {
  return(length(vect[which(vect == text)]))
}
lift = lift[-which(lift$V4 %in% traits_to_remove),]
lift = lift[which(lift$V4 %in% traits_to_keep),]

#lift = lift[which(lift$V4 == "Alzheimer.s.disease"),]

getOverlapP <- function(peaks, lift) {
  all_traits = unique(lift$V4)
  #coors = rownames(peaks)
  chr = peaks$V1
  start = peaks$V2 
  end = peaks$V3 
  gr1 <- GRanges(seqnames=chr,
                 ranges=IRanges(start,end),
                 strand = "+")
  snps = GRanges(seqnames=lift$V1,
                 ranges=IRanges(lift$V2,lift$V3),
                 strand = "+", snp = lift$V5, info = lift$V4, tag = lift$V6)
  ov = findOverlaps(gr1, snps,
                    maxgap=-1L, minoverlap=0L,
                    type=c("any", "start", "end", "within", "equal"),
                    select=c("all", "first", "last", "arbitrary"),
                    ignore.strand=T)
  overlap_table = data.frame(gr1[ov@from])
  overlap_table = cbind(overlap_table, data.frame(snps[ov@to]))

  over_count = c()
  for( t in all_traits) {
    # over_count= length(overlap_table$info[which(overlap_table$info == t)])
    over_count = c(over_count, countOccurance(t , overlap_table$info))
  }
  names(over_count) = all_traits
  #snp_count = table(lift$V4)[all_traits]
  background = data.frame(fread("~/projects/mouse_aging/analysis/snapATAC/all_celltypes/all_celltypes.cluster_peaks.merged.bed"))
  gr2 = GRanges(seqnames=background$V1,
               ranges=IRanges(background$V2,background$V3),
               strand = "+")
  ov2 = findOverlaps(gr2, snps,
                    maxgap=-1L, minoverlap=0L,
                    type=c("any", "start", "end", "within", "equal"),
                    select=c("all", "first", "last", "arbitrary"),
                    ignore.strand=T)
  overlap_table2 = data.frame(gr2[ov2@from])
  overlap_table2 = cbind(overlap_table2, data.frame(snps[ov2@to]))
  snp_count = c()
  for( t in all_traits) {
    # over_count= length(overlap_table$info[which(overlap_table$info == t)])
    snp_count = c(snp_count, countOccurance(t , overlap_table2$info))
  }
  names(snp_count) = all_traits
  
  test_width = sum(peaks$V3-peaks$V2)
  background_width = sum(background$V3-background$V2)
  odds = c()
  pvals = c()
  for(t in all_traits) {
    if (snp_count[t] > 0 ){ 
    test = binom.test(x=as.numeric(over_count[t]), n=snp_count[t], p =(test_width/background_width),
                      alternative = c("greater"),
                      conf.level = 0.95)
    odds = c(odds, over_count[t]/snp_count[t]/ (test_width/background_width))
    pvals = c(pvals, test$p.value)
    }
    else {
      odds = c(odds, 1)
      pvals = c(pvals,1)
    }
  }
 # print(over_count,pvals)
  return(list(over_count,odds, pvals))
}


all_peak_matrix = data.frame(unique(lift$V4))
counts = data.frame(unique(lift$V4))
odds = data.frame(unique(lift$V4))

for( i in 1:nrow(key)) {
 tissue = key$Tissue[i]
 cl = key$Cluster[i]
 print(c(tissue,cl))
 peaks = data.frame(fread(paste("/mnt/tscc/yanxiao/projects/mouse_aging/analysis/snapATAC/", tissue,"/peak.cluster/",tissue,".metacell_",cl,"_peaks.narrowPeak", sep = ""), header = F, ))
 #rownames(peaks) = paste0(peaks$V1,":",peaks$V2,"-",peaks$V3)
 peaks = peaks[,c(1:3)]
 
 #peaks = peaks[order(rowMeans(peaks[,1:6]), decreasing = T),]
 #peaks = peaks[order(peaks$logCPM,decreasing=T),]
# write.table(sub("(.*):(.*)-(.*)","\\1\t\\2\t\\3",rownames(peaks)), 
#             paste0("~/projects/mouse_aging/tests/GWAS/",tissue,".",key[i,"Name"],".peaks.bed"),row.names=F,quote=F,col.names=F)
 
 out= getOverlapP(peaks = peaks,lift = lift)
 counts = cbind(counts, out[[1]])
 odds = cbind(odds,out[[2]])
 all_peak_matrix = cbind(all_peak_matrix, out[[3]])
   
 }

colnames(counts)[-1] = colnames(odds)[-1] =colnames(all_peak_matrix)[-1] = paste(key[,"Tissue",],key[,"Name"],sep=".")
rownames(all_peak_matrix) = all_peak_matrix[,1]
all_peak_matrix = all_peak_matrix[,-1]
rownames(odds) = odds[,1]
odds = odds[,-1]
all_tog_p = all_peak_matrix

annotation_col_all_tog = data.frame(
  Tissue = key$Tissue,
  Clade = key$Clade
)

#got the cell type ordering. 
#cts = read.table("/mnt/tscc/yanxiao/projects/mouse_aging/analysis/snapATAC/all_celltypes/Spearman.clustering.celltypes.ordered.txt",stringsAsFactors=F)$V1
#cts = sub("(..)\\...?\\.(.*)","\\1.\\2",cts)
#all_tog_p = all_tog_p[,match(cts,colnames(all_tog_p))]

rownames(annotation_col_all_tog) = colnames(all_tog_p)
pdf("/mnt/tscc/lamaral/projects/Aging/aging_share/figures/GWAS/enrichment.pdf",height = 10, width = 15)

pheatmap(-log10(all_tog_p[which(rowMins(data.matrix(all_tog_p))<0.01),!is.na(key$Clade)]),display_numbers = T, #cluster_cols =F,
#  clustering_distance_rows = "correlation",
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), breaks = seq(1,8,length.out = 100), number_format  = "%.0f",
         annotation_col = annotation_col_all_tog,main = "All together GWAS enrichment" )
dev.off()

#pheatmap(odds,color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), breaks = seq(1,3,length.out = 100), number_format  = "%.0f",)
#tab1 = t(counts[-1])
#write.table(tab1,"~/projects/mouse_aging/tests/GWAS/overlap_cnts.txt")
#tab1 = t(all_peak_matrix[-1])
#write.table(tab1,"~/projects/mouse_aging/tests/GWAS/overlap_pval.txt")


