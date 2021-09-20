# colnames(a2) = c("H3K27ac","H3K27me3","H3K4me3","H3K4me1","H3K9me3","CTCF")

setwd("../../analysis/all_celltypes.all_tissue")
a = data.frame(fread("allcelltypes.cluster_peaks.merged.thresh.annotated.txt"))
as = colSums(a[-c(1:3)])
b = data.frame(fread("control.merged.thresh.annotated.txt"))
bs = colSums(b[-c(1:3)])

dat = data.frame(marks=c("H3K27ac","H3K27me3","H3K4me3","H3K4me1","H3K9me3","CTCF"),fold=as/bs)

pdf("Enrichment_overlap_with_marks.thresh.pdf")
ggplot(dat) + geom_col(aes(x=marks, y=fold)) +
  geom_hline(yintercept=1,color="black",linetype="dashed") + 
  coord_flip() +
  theme_classic(base_size=20)
dev.off()

