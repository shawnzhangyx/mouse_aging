################## count the number of differential peaks using the same cut-off.

setwd("../../analysis/snapATAC/DH/age_diff_edgeR.snap")
## read the meta info. 
meta = read.delim("../DH.pool.barcode.meta_info.txt")



files = list.files(pattern=".edger.txt")
max_cluster = length(files)
dat = list()


for (cl in 1:max_cluster) {
dat[[cl]] = read.csv(paste0(cl,".edger.txt"))
dat[[cl]]$cluster = cl
}

dat2 = do.call(rbind, dat)

t2 = table(dat2$PValue < 0.001, dat2$logFC>0, dat2$cluster)
mt2 =melt(t2)

ncell = table(meta$cluster)
mt2$ncells = ncell[match(mt2$Var3,names(ncell))]

pdf("number_of_diff_peaks.pdf")
ggplot(subset(mt2,Var1==TRUE)) + 
  geom_col(aes(x=factor(Var3,levels=(max_cluster:1)),
      y=value,fill=Var2),position="dodge") + 
  scale_fill_discrete(name="Up in Aging") + 
  coord_flip() 

ggplot(subset(mt2,Var1==TRUE)) + geom_text(aes(y=value,x=ncells,label=Var3,color=Var2)) +
  scale_color_discrete(name="Up in Aging") + 
  xlab("Number of cells") +
  ylab("Number of Differential Peaks") + 
  coord_flip()


dev.off()


# Number of Shared differential peaks.
diff = dat2[dat2$PValue < 0.001,]
diff_peak_cnt = sort(table(diff$X))
shared_diff = data.frame(table(table(diff$X))[-1])

pdf("Diff_peaks_shared_by_clusters.pdf")
ggplot(shared_diff) + geom_col(aes(x=Var1,y=Freq)) +
  geom_text(aes(x=Var1,y=Freq,label=Freq),nudge_y=100)+
  xlab("Number of Cell Types found")
dev.off()

mat = matrix(0,nrow=max_cluster,ncol=max_cluster)

for (i in 1:max_cluster){
  for (j in 1:max_cluster){
    print(paste(i,j))
    set1 = diff[diff$cluster==i,"X"] 
    set2 = diff[diff$cluster==j,"X"]
    mat[i,j] = length(intersect(set1,set2))/length(union(set1,set2))
    }}

melted = melt(mat)

pdf("overlap_of_differential_peaks.between_cluster.pdf")
ggplot(subset(melted,Var1!=Var2)) + 
  geom_tile(aes(x=Var1,y=Var2,fill=value)) + 
  scale_fill_gradient(high="red",low="white",name="Jaccard Score")
dev.off()





