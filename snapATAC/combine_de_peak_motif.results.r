tissue=commandArgs(trailing=T)[1]

setwd(paste0("../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/motif.homer"))

files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

dat = list()
for (file in files){
  cluster = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\1",file)
  category = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\2",file)
  tmp = read.delim(file)
  dat[[file]] = data.frame(cluster=cluster,category=category,tmp[1:5,c(1,2,3)])
  }

out = do.call(rbind,dat)
write.csv(subset(out,category!="both"),"DE_Peaks.no_background.homer_results.csv")

#setwd("../../analysis/snapATAC/DH/age_diff_edgeR.snap/motif.homer.bg")
setwd("../motif.homer.bg/")

files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

dat = list()
for (file in files){
  cluster = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\1",file)
  category = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\2",file)
  tmp = read.delim(file)
  dat[[file]] = data.frame(cluster=cluster,category=category,tmp[1:5,c(1,2,3)])
  }

out = do.call(rbind,dat)
write.csv(subset(out,category!="both"),"DE_Peaks.with_background.homer_results.csv")

setwd("../motif.homer.csbg/")

files = list.files(pattern = "knownResults.txt",full.names=T,recursive=T)

dat = list()
for (file in files){
  cluster = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\1",file)
  category = sub("\\.\\/(.*)\\.(.*)\\.homer.*","\\2",file)
  tmp = read.delim(file)
  dat[[file]] = data.frame(cluster=cluster,category=category,tmp[1:5,c(1,2,3)])
  }

out = do.call(rbind,dat)
write.csv(subset(out,category!="both"),"DE_Peaks.with_csbackground.homer_results.csv")

