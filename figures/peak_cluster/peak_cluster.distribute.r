path = "/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/clustered_diff_peak_locations"

files = list.files(pattern=".*clustered.*.txt",path,recursive=T)
files = files[!grepl("old_May5",files)]
files = files[!grepl("old_2020May6",files)]
files = files[grepl("(DH|FC)",files)]


dat = list()

for (f in files){
tmp = read.delim(paste0(path,"/",f))
if (nrow(tmp) > 0 ){
 dat[[f]] = tmp
 dat[[f]]$tissue = sub("(..)\\/.*peaks(.*)_0.001.txt","\\1",f)
 dat[[f]]$direct = sub("(..)\\/.*peaks(.*)_0.001.txt","\\2",f)
 }
}

dat2 = do.call(rbind,dat)
dat2$start = as.numeric(sub(".*:(.*)-(.*)","\\1",dat2$locations))
dat2$end = as.numeric(sub(".*:(.*)-(.*)","\\2",dat2$locations))
dat2$size = dat2$end-dat2$start


ggplot(dat2) + geom_point(aes(x=size,y=counts,color=smth20_score,size=smth20_score))

dat2 = dat2[order(-dat2$counts),]
dat2 = dat2[order(-dat2$size),]
dat2 = dat2[order(-dat2$smth20_score),]


ggplot(dat2) + geom_point(aes(x=counts,y=smth20_score))


dat3= dat2[which(dat2$smth20_score>0.4 & dat2$size>100000),]
write.table(dat3, "clustered_diff_peak_all.smth20_gt0.4.size_gt100k.txt",quote=F,sep="\t",row.names=F)

dat3$chr = sub("(.*):.*","\\1",dat3$locations)

out = dat3[,c("chr","start","end")]
out$TisCluDir = paste0(dat3$tissue,".",dat3$clust,".",dat3$direct)
out = out[order(out$chr,out$start,out$end),]
write.table(out,"clustered_diff_peak.bed",row.names=F,col.names=F,sep="\t",quote=F)

