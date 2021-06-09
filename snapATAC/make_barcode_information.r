tissue= commandArgs(trailing=T)[1]
input = commandArgs(trailing=T)[2]

setwd(paste0("../../analysis/snapATAC/",tissue))

#a=read.delim(paste0("snapFiles/",tissue,".pool.snapATAC.TSS_filter.cluster.meta.txt"))
a=read.delim(paste0("snapFiles/",input))

# read barcode information
bfiles = list.files("../../../data/snATAC/qc/bc_info_by_sample/",pattern=tissue)

ll = list()
for (idx in 1:6){
tmp=read.delim(paste0("../../../data/snATAC/qc/bc_info_by_sample/",bfiles[idx]))
tmp$sample=substr(bfiles[idx],1,10)
ll[[idx]] = tmp
}

b = do.call(rbind,ll)

a$TN = b$raw[match(paste(a$sample,a$barcode),paste(b$sample,b$barcode))]
a$MT.pc = b$mt.pc[match(paste(a$sample,a$barcode),paste(b$sample,b$barcode))]
#a$FIP.pc = b$FIP.pc[match(paste(a$sample,a$barcode),paste(b$sample,b$barcode))]
#a$TSS.pc = b$TSS.pc[match(paste(a$sample,a$barcode),paste(b$sample,b$barcode))]

# some adjustment of column names. 
a$sample=NULL
#colnames(a)[match(c("x.sp.sample","x.sp.cluster"),colnames(a))]= c("sample","cluster")
# after harmony
colnames(a)[match(c("x.after.sp.sample","x.after.sp.cluster"),colnames(a))]= c("sample","cluster")
a$stage = substr(a$sample,4,5)
a$rep = substr(a$sample,7,10)
a$landmark = NULL

write.table(a,paste0(tissue,".pool.barcode.meta_info.txt"),row.names=F,quote=F,
  sep="\t")



