sample=commandArgs(trailing=T)[1]
setwd("../../data/snATAC/qc")

raw = read.table(paste0("../bam.raw/",sample,"/",sample,".raw.barcode.cnts.txt"))
dedup = read.table(paste0("../bam.dedup/",sample,"/",sample,".dedup.barcode.cnts.txt"))
filter = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.cnts.txt"))
mt = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.chrM.cnts.txt"))
FIP = read.table(paste0("../ct2peaks/",sample,".FIP.barcodes.cnts.txt"))


out = Reduce(function(...)merge(...,by="V1"),list(raw,dedup,filter,mt,FIP))

colnames(out) = c("barcodes","raw","dedup","filter","mt","FIP")
out = out[order(-out$filter),]
out$mt.pc = out$mt/out$filter*100
out$FIP.pc = out$FIP/out$filter*100


write.table(out,paste0("bc_info_by_sample/",sample,".barcode.info.txt"),row.names=F,sep="\t",quote=F)

