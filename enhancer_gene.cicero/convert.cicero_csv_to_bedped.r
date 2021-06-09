args = commandArgs(trailing=T)
input = args[1]
output = args[2]

a=read.csv(input) 
peak1 = gsub("_","\t",a$Peak1)
peak2 = gsub("_","\t",a$Peak2)

write.table(cbind(peak1,peak2),output,row.names=F,col.names=F,quote=F,sep="\t")
