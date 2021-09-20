setwd("../../analysis/all_celltypes.all_tissue/")

a=read.table"cluster.read_counts.txt")


a$V3 = ifelse(a$V2<5e6,     0.1, 
          ifelse(a$V2<25e6, 0.05,
          ifelse(a$V2<50e6, 0.025,
          ifelse(a$V2<100e6,0.01,
                            0.001))))

write.table(a,"cluster.read_counts.fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")

