#Endo:
# DH.13
# FC.15
# LM.6
# HT.7
# HT.3
# BM.15

#all 27. 

a=read.delim("/projects/ps-renlab/lamaral/projects/Aging/all_tissues/all_merged_40kL.cluster.meta.txt")

#files = list.files(path="/mnt/silencer2/home/shz254/projects/mouse_aging/analysis/snapATAC/",pattern=".meta_info.txt",recursive=T,full.name=T)


#dat_list= list()
#for (file in files) {
#  b = read.delim(file)
#  dat_list[[file]] = b
#  }

#dat = do.call(rbind,dat_list)
#row.names(dat) = NULL

#b = dat[,c("sample","barcode","cluster")]

#d = merge(a,b,by=c("sample","barcode"))
#d$tissue_cluster= paste0(d$tissue,"_",d$cluster_y)

#table(d$cluster.x,d$tissue_cluster)

#table(a$cluster,a$tissue_cluster)[27,]

a2= a[which(a$tissue_cluster %in% c("DH_13","FC_15","LM_6","HT_3","HT_7","BM_15")),]
a2$tissue_cluster = as.character(a2$tissue_cluster)

