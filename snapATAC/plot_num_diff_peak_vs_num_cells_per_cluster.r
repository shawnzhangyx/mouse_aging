##
library(ggrepel)
#tissue=commandArgs(trailing=T)[1]

setwd(paste0("../../analysis/snapATAC/"))

out = list()

for ( tissue in c("HT","DH","FC","LM")) {
diff = read.delim(paste0(tissue,"/age_diff_bycluster/",tissue,".peaks.age_diff.summary.txt"))
meta = read.delim(paste0(tissue,"/snapFiles/",tissue,".pool.snapATAC.cluster.meta.txt"))

meta$cluster=paste0("C",meta$x.sp.cluster)

num_cells1 = table(meta$cluster)
num_reads = aggregate(UM~cluster,data=meta,sum)


diff$num_cells = num_cells1[match(diff$cluster,names(num_cells1))]
diff = merge(diff,num_reads,by="cluster")
diff$cluster= paste0(tissue,".",diff$cluster)
diff$tissue = tissue

out[[length(out)+1]] = diff

pdf(paste0(tissue,"/",tissue,".num_diff_peak.vs.num_cells_reads.pdf"))

print( ggplot(diff) + geom_point(aes(x=Down_Reg+Up_Reg,y=num_cells)) + 
  geom_text_repel(aes(x=Down_Reg+Up_Reg,y=num_cells,label=cluster)))

print( ggplot(diff) + geom_point(aes(x=Down_Reg+Up_Reg,y=UM)) +
  geom_text_repel(aes(x=Down_Reg+Up_Reg,y=UM,label=cluster)) + 
  scale_y_log10())
dev.off()

}

dat = do.call(rbind,out)
pdf("all_tissues..num_diff_peak.vs.num_cells_reads.pdf")
ggplot(dat) + geom_point(aes(x=Down_Reg+Up_Reg,y=num_cells)) +
  geom_text_repel(aes(x=Down_Reg+Up_Reg,y=num_cells,label=cluster,color=tissue)) 

ggplot(dat) + geom_point(aes(x=Down_Reg+Up_Reg,y=UM)) +
  geom_text_repel(aes(x=Down_Reg+Up_Reg,y=UM,label=cluster,color=tissue)) +
    scale_y_log10()
dev.off()
#ggplot(diff) + geom_point(aes(x=Down_Reg,y=Up_Reg))

