

setwd("../../data/snATAC/")

sample_info = read.table("../../scripts/pre_processing/demultiplex.sample_info.txt")


#total_reads = NULL
bamRaw = NULL
bamDedup = NULL
fracDup = NULL
bamFilter = NULL

samples=rownames(sample_info)
for (sample in samples){
  print(sample)
#  fastq = fread(paste0("fastq.demultiplex/",sample,"/", sample ,".barcode.cnts.txt"))
#  total_reads = c(total_reads, sum(fastq$V2))
  bam = fread(paste0("bam.raw/",sample,"/" , sample , ".raw.barcode.cnts.txt"))
  bamRaw = c(bamRaw, sum(bam$V2))
  # the 8th line of dedup QC file.   
  dedup = readLines(paste0("bam.dedup/",sample,"/" , sample , ".dedup.qc"))[8]
  nums = as.numeric(strsplit(dedup,split="\t")[[1]][-1])
  bamDedup = c(bamDedup,nums[1]+nums[2]-nums[5]-nums[6])
  fracDup = c(fracDup,nums[8])
#  bamDedup = c(bamDedup, sum(bam$V2))
  bam = fread(paste0("bam.filter/",sample,"/" , sample , ".filter.barcode.cnts.txt"))
  bamFilter = c(bamFilter, sum(bam$V2))
}


out = data.frame(samples,bamRaw,bamDedup,fracDup,bamFilter)
out$fracFilter = 1-out$bamFilter/out$bamDedup

write.table(out,"qc/Sample_ReadTotal.txt",quote=F,row.names=F,sep="\t")
#out$Dup = 


