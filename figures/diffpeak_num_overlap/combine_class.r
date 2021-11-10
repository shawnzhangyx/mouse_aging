a=read.delim("../celltype_annotation.txt")

clades = unique(a$Clade)
clades= clades[!is.na(clades)]
for (clade in clades ) {
  tmp = a[which(a$Clade==clade),]
  dat.list = list()
  for (i in 1:nrow(tmp)){
    tissue=tmp[i,"Tissue"]
    cluster=tmp[i,"Cluster"]
    fname = paste0("../../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/",cluster,".up.bed")
    dat.list[[fname]] = fread(fname)
    fname = paste0("../../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/",cluster,".down.bed")
    dat.list[[fname]] = fread(fname)
    }
  dat = do.call(rbind,dat.list)
  write.table(dat, paste0("class_diff_peak/",clade,".bed"),row.names=F,col.names=F,sep="\t",quote=F)
    }


