tissue="DH"


for (tissue in c("DH","FC","HT","LM","BM")) {
  print(tissue)
  a=read.delim(paste0("../../../analysis/snapATAC/",tissue,"/",tissue,".pool.barcode.meta_info.txt"))
a$cluster = factor(a$cluster)

print(paste0(tissue,".umap.cluster.png"))
png(paste0(tissue,".umap.cluster.png"),height=1080,width=1080,res=300)
print(ggplot(a) + geom_point(aes(x=umap.1,y=umap.2,color=cluster),size=0.1,alpha=0.1)+
  theme_classic() +
  theme(legend.position="none") )

  dev.off()

}
