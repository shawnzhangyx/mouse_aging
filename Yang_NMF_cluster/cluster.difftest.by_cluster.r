tissue=commandArgs(trailing=T)[1]
rank=commandArgs(trailing=T)[2]
setwd(paste0("../../analysis/Yang_NMF_method/",tissue,"/R",rank))
system("mkdir age_diff_bycluster_specific_peak")

library(edgeR)
library(doParallel)
registerDoParallel(cores=rank)

summ = foreach(grp=1:rank) %dopar% {

a=read.delim(paste0("counts/",tissue,".C",grp,".counts"),skip=1)
b =read.delim(paste0("counts/",tissue,".C",grp,".counts.summary"))
cnts = a[,-c(1:6)]
total = colSums(b[c(1,10),match(colnames(cnts),colnames(b))])

samples = paste0("C",sub(".*metacell_(.{1,2}\\....rep.).bam","\\1",colnames(cnts)))
#grps = paste0("C",sub(".*metacell_(.{1,2})\\....rep..bam","\\1",colnames(cnts)))
stage = sub(".*metacell_.{1,2}\\.(..).rep..bam","\\1",colnames(cnts))

colnames(cnts) = samples
rownames(cnts) = paste0(a$Geneid,":",a$Chr,":",a$Start,"-",a$End)

y = DGEList(cnts)
#y$sample$lib.size=total[idx]
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)

print(grp)
groups = stage
design = model.matrix(~0+groups)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

contrast.matrix = matrix(c(c(-1,1,0),c(0,-1,1),c(-1,0,1)),nrow=3)
lrt = glmLRT(fit_tag, contrast =contrast.matrix)
fdr = p.adjust(lrt$table$PValue,method="BH")
out = cbind(cpm(y),lrt$table,fdr)
out = out[order(out$PValue),]


write.table(out,paste0("age_diff_bycluster_specific_peak/",tissue,".",grp,".peaks.age_diff.txt"),sep="\t",quote=F)

sig = factor(out$logFC.3[out$fdr<0.05]>0,levels=c("FALSE","TRUE"))

tab = table(sig)
if (sum(tab) ==0) { tab = c(0,0) }

up = rownames(out)[which(out$logFC.3>0 & out$fdr<0.05)]
up = sub("peak.*:(.*):(.*)-(.*)","\\1\t\\2\t\\3",up)
 if (length(up) > 50) { write(up, paste0("age_diff_bycluster_specific_peak/",tissue,".",grp,".peaks.up.bed")) }
down = rownames(out)[which(out$logFC.3<0 & out$fdr<0.05)]
down = sub("peak.*:(.*):(.*)-(.*)","\\1\t\\2\t\\3",down)
 if (length(down) > 50) { write(down, paste0("age_diff_bycluster_specific_peak/",tissue,".",grp,".peaks.down.bed")) }
comb = c(up,down) 
 if (length(comb) > 50) { write(comb, paste0("age_diff_bycluster_specific_peak/",tissue,".",grp,".peaks.comb.bed")) }


c(grp,tab)
}

summ2 = do.call(rbind,summ)
colnames(summ2) = c("cluster","Down_Reg","Up_Reg")
summ2 = summ2[order(as.numeric(sub("C(.*)","\\1",summ2[,1]))),]
write.table(summ2,paste0("age_diff_bycluster_specific_peak/",tissue,".peaks.age_diff.summary.txt"),sep="\t",quote=F,row.names=F)
