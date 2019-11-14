sample=commandArgs(trailing=T)[1]
setwd("../../data/snATAC/qc")

## raw,mapped, deduplicated, filtered, mitochondria reads. 
raw = read.table(paste0("../fastq.demultiplex/",sample,"/",sample,".barcode.cnts.txt"))
#mapped = read.table(paste0("../bam.raw/",sample,"/",sample,".raw.barcode.cnts.txt"))
dedup = read.table(paste0("../bam.dedup/",sample,"/",sample,".dedup.barcode.cnts.txt"))
filter = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.cnts.txt"))
mt = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.chrM.cnts.txt"))
# Fragments in TSS. 
tss_read = read.table(paste0("../tss_barcode/",sample,".barcode.TSS.cnts.txt"))
# Fragments in peak. 
FIP = read.table(paste0("../ct2peaks/",sample,".FIP.barcodes.cnts.txt"))

# TSS enrichment score. 
tss_enrich = read.table(paste0("/projects/ps-renlab/lamaral/projects/Aging/TSS_enrich/",sample,"_TSS_enrich_res.txt"))[,c(1,4)]

out = Reduce(function(...)merge(...,by="V1"),list(raw,dedup,filter,FIP,tss_enrich))
out$mt = mt$V2[match(out$V1,mt$V1)]
out$mt[is.na(out$mt)] = 0

out$tss_read = tss_read$V2[match(out$V1,tss_read$V1)]
out$tss_read[is.na(out$tss_read)] = 0


colnames(out) = c("barcodes","raw","dedup","filter","FIP","TSS_enrich","mt","TSS_read")
out = out[order(-out$filter),]
out$mt.pc = out$mt/out$filter*100
out$FIP.pc = out$FIP/out$filter*100
out$TSS.pc = out$TSS_read/out$filter*100

write.table(out,paste0("bc_info_by_sample/",sample,".barcode.info.txt"),row.names=F,sep="\t",quote=F)

