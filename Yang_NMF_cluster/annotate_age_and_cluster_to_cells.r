#tissue=commandArgs(trailing=T)[1]
statH = commandArgs(trailing=T)[1]
outFile=commandArgs(trailing=T)[2]

#setwd(paste0("../../analysis/Yang_NMF_method/",tissue))
a=read.table(statH)


a$V3 = a$V3+1
a$stage = substr(a$V1,1,2)
a$rep = substr(a$V1,4,7)
a$barcode=substr(a$V1,9,100)

out = a[,c(10,3,8,9)]
colnames(out) = c("barcode","cluster","stage","rep")
write.table(out,outFile,row.names=F,sep="\t",quote=F)

