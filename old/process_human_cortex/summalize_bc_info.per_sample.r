sample=commandArgs(trailing=T)[1]
setwd("../../data/human_cortex/qc")

## raw,mapped, deduplicated, filtered, mitochondria reads. 
raw = read.table(paste0("../fastq.relabel/",sample,".barcode.cnts.txt"))
#mapped = read.table(paste0("../bam.raw/",sample,"/",sample,".raw.barcode.cnts.txt"))
dedup = read.table(paste0("../bam.dedup/",sample,"/",sample,".dedup.barcode.cnts.txt"))
filter = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.cnts.txt"))
mt = read.table(paste0("../bam.filter/",sample,"/",sample,".filter.barcode.chrM.cnts.txt"))
# Fragments in TSS. 

# TSS enrichment score. 
tss_enrich = read.table(paste0("../TSSvsFrag.out/",sample,"/stat.txt"))[,c(1,4)]

out = Reduce(function(...)merge(...,by="V1"),list(raw,dedup,filter,tss_enrich))
out$mt = mt$V2[match(out$V1,mt$V1)]
out$mt[is.na(out$mt)] = 0



colnames(out) = c("barcodes","raw","dedup","filter","TSS_enrich","mt")
out = out[order(-out$filter),]
out$mt.pc = out$mt/out$filter*100
#out$TSS.pc = out$TSS_read/out$filter*100

write.table(out,paste0("bc_info_by_sample/",sample,".barcode.info.txt"),row.names=F,sep="\t",quote=F)

