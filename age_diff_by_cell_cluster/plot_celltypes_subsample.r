setwd("../../analysis/snapATAC/de_peaks.subsample")
meta = read.delim("../../../aging_share/figures/celltype_annotation.txt")
reads= read.delim("all_celltypes.all_sample.counts.summary")
reads = data.frame(colSums(reads[,-1]))
colnames(reads) = "reads"
reads$sample = sub(".*\\.(..)\\.metacell_(.*)\\.sorted.bam","\\1.\\2",rownames(reads))
reads$ct = sub("(..)\\.(.*)\\.(..)\\.(rep.)","\\1.\\2",reads$sample)
reads$age= sub("(..)\\.(.*)\\.(..)\\.(rep.)","\\3",reads$sample)
reads$class = meta$Clade[match(reads$ct, paste0(meta$Tissue,".",meta$Cluster))]

reads$name = paste0(meta$Tissue,".",meta$Name)[match(reads$ct, paste0(meta$Tissue,".",meta$Cluster))]
#ggplot(reads) + geom_boxplot(aes(x=class,y=log10(reads)))
filt = reads[which(reads$reads>5e5 & reads$age !="10"),]
cts = names(table(filt$name)[table(filt$name)==4])
write(cts)
filt = reads[which(reads$reads>1e6 & reads$age !="10"),]
cts = names(table(filt$name)[table(filt$name)==4])
write(cts)
filt = reads[which(reads$reads>2e6 & reads$age !="10"),]
cts = names(table(filt$name)[table(filt$name)==4])
write(cts)

filt = reads[which(reads$reads>5e6 & reads$age !="10"),]
cts = names(table(filt$name)[table(filt$name)==4])
write(cts)

