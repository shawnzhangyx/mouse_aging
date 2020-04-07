a=data.frame(fread("/projects/ps-renlab/yanxiao/software/RepEnrich2/refs/setup_folder_mm10/repnames.bed"))

b=read.delim("../../analysis/repeat_analysis/repeat_class_family.wo_simple_repeat.txt",header=F)

a = a[which(a$V4 %in% b$V1),]
a$len = a$V3-a$V2

num = data.frame(table(a$V4))
num = num[order(-num$Freq),]
num$class = b$V2[match(num$Var1,b$V1)]

len = aggregate(len~V4,a,sum)
len = len[order(-len$len),]
len$class = b$V2[match(len$V4,b$V1)]

write.table(num,"../../analysis/repeat_analysis/summary/repeat_num.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(len,"../../analysis/repeat_analysis/summary/repeat_cov.txt",row.names=F,col.names=F,sep="\t",quote=F)

