#celltypes = read.csv("All_tissues.diffpeak_share.table.p0.01.csv")
celltypes=read.table("celltypes.ordered_heatmap.txt",stringsAsFactors=F)$V1
#celltypes = colnames(celltypes)[-1]

meta = read.delim("../celltype_annotation.txt")

ct = paste0(meta$Tissue,".",meta$Cluster)[match(celltypes,paste0(meta$Tissue,"_",meta$Name))]


dat_list = list()
for (  tissue in c("DH","FC","HT","LM","BM") ){
  fs = list.files(pattern="edger.txt",path=paste0("../../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap"),full.names=T)
  for ( n in 1:length(fs)) {
    diff = data.frame(fread(fs[n]))
    cl = as.numeric(sub(".*snap\\/(.*).edger.txt","\\1",fs[n]))
    len = length(which(diff$PValue < 0.01))
    dat_list[[length(dat_list)+1]] = data.frame(tissue,cl,len)
}
}
diff = do.call(rbind,dat_list)
diffnum = diff$len[match(ct, paste0(diff$tissue,".",diff$cl))]

a = data.frame(od=celltypes,diffnum=diffnum)
a$od = factor(a$od,levels=rev(celltypes))

g3 =
  ggplot(a) + geom_col(aes(x=od,y=diffnum)) +
  geom_text(aes(x=od,y=diffnum,label=diffnum),hjust=-0.2) +
  coord_flip() + ylim(0,5000) +
  theme_bw()

pdf("selected_celltypes.diffnum.pdf",width=4,height=15)
g3
dev.off()

