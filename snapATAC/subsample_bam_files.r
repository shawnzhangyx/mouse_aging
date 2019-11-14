setwd("../../analysis/snapATAC/DH/"

files = list.files(pattern="sorted.bam$",path="bam.cluster",full.names=T)

file = files[1]

library(doParallel)
registerDoParallel(cores=length(files))

out = foreach(file=files)%dopar% {
system(paste("samtools view -c", file),intern=T)
}
size = as.numeric(out)
frac = min(size)/size

system("mkdir bam.cluster.subsample")

foreach( idx=1:length(size) )%dopar% {
 if ( frac[idx] == 1 ){
    system(paste0("samtools view -f 64 ", files[idx], sub("bam.cluster"," -o  bam.cluster.subsample",files[idx])))
    } else {
    system(paste("samtools view -f 64 -s ", frac[idx]+1, files[idx], sub("bam.cluster","-o bam.cluster.subsample",files[idx]) ) )
  }

}


### call peaks 

foreach( idx=1:length(size) )%dopar% {
  system(paste("macs2 callpeak -t ", sub("bam.cluster","bam.cluster.subsample",files[idx]), sub("bam.cluster"," -n bam.cluster.subsample",files[idx]))) 
  }


peaks = list.files(pattern="narrowPeak",path="bam.cluster.subsample",full.names=T)
out2 = foreach(file=peaks)%dopar% {
system(paste("wc -l", file),intern=T)
}
out2 = do.call(rbind,out2)
numP = as.numeric(sub("(.*) .*","\\1",out2))
clu = as.numeric(sub(".*metacell_(.*).sorted.*","\\1",peaks))
tab = data.frame(clu,numP)
pdf("bam.cluster.subsample/num_peaks_per_cluster.pdf")
ggplot(tab) + geom_col(aes(factor(clu),numP))
dev.off()



