#tissue=commandArgs(trailing=T)[1]
statH = commandArgs(trailing=T)[1]
outFile=commandArgs(trailing=T)[2]

#setwd(paste0("../../analysis/Yang_NMF_method/",tissue))
a=read.table(statH,header=F,skip=1)


#a$V3 = a$V3+1
a$stage = substr(a$V1,4,5)
a$rep = substr(a$V1,7,10)
#a$barcode=a$V2#substr(a$V1,12,100)
# cluster V8

out = a[,c(2,8,11,12)]
colnames(out) = c("barcode","cluster","stage","rep")
write.table(out,outFile,row.names=F,sep="\t",quote=F)

