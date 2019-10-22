a=fread("/projects/ps-renlab/yanxiao/software/RepEnrich2/refs/mm10_repeatmasker.txt")
b=a[,c(10,11)]
b=b[-which(duplicated(b$V10)),]
b=as.data.frame(b)
b$family = sub("(.*)\\/.*","\\1",b$V11)
write.table(b,"../../analysis/repeat_analysis/repeat_class_family.txt",quote=F,sep="\t",row.names=F,col.names=F)
b = b[!b$family=="Simple_repeat",]
write.table(b,"../../analysis/repeat_analysis/repeat_class_family.wo_simple_repeat.txt",quote=F,sep="\t",row.names=F,col.names=F)

rep1 = fread("../../analysis/repeat_analysis/repnames.txt",header=F)
rep1 = rep1[substr(rep1$V1,1,1)!="_",]
write.table(rep1,"../../analysis/repeat_analysis/repnames.wo_simple_repeat.txt",row.names=F,col.names=F,quote=F)

