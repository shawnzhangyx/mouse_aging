tissue="DH"
setwd(paste0("../../analysis/snapATAC/",tissue))

a=read.delim("DH.peaks.counts",skip=1)
cnts = a[,-c(1:6)]
#total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])
total = colSums(cnts)
med = median(total)
cnts = sweep(cnts,2,total,'/')*max(total)

samples = paste0("C",sub(".*metacell_(.{1,2}\\....rep.).sorted.bam","\\1",colnames(cnts)))
grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..sorted.bam","\\1",colnames(cnts)))
grps = factor(grps,levels=unique(grps))
stage = sub(".*metacell_.{1,2}\\.(..).rep..sorted.bam","\\1",colnames(cnts))


reads = as.numeric(cnts[grep("peak_112049",a$Geneid),])
dat = data.frame(grps,reads)

ggplot(dat) + geom_boxplot(aes(grps,reads))

