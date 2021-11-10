mat = read.csv("All_tissues.diffpeak_share.table.p0.01.csv",row.names=1)
rownames(mat) = colnames(mat) = sub("_",".",colnames(mat))


sim = read.csv("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/snapATAC/all_celltypes/Spearman.correlation.celltypes.csv",row.names=1)
rownames(sim) = colnames(sim) = sub("(..)\\.(.*?)\\.(.*)","\\1.\\3",colnames(sim))


#sim.melt = melt(as.matrix(sim))
sim = sim[which(rownames(sim) %in% rownames(mat)), which( colnames(sim) %in% rownames(mat))]


dat_list = list()
for ( x in rownames(mat)) {
  p = sort(sim[which(rownames(sim)==x),],decreasing=T)
  tmp = data.frame(ct1 = x,ct2=names(p),dist=as.numeric(p) ,distR = 1:length(p))
  dat_list[[x]] = tmp

 }

dat = do.call(rbind,dat_list)

mat.melt = melt(as.matrix(mat))
colnames(mat.melt) = c("ct1","ct2","jac")

dat2 = merge(dat,mat.melt,by = c("ct1","ct2"))

ggplot(subset(dat2,distR>1)) + geom_boxplot(aes(x=factor(distR),y=jac))

#ggplot(subset(dat2,distR>1)) + geom_violin(aes(x=factor(distR),y=jac))

pdf("jaccard_distance_by_celltype_similarity.pdf")
ggplot(subset(dat2,distR>1)) + 
  geom_point(aes(x=1-dist,y=jac),color="grey") +
  geom_smooth(aes(x=1-dist,y=jac)) +
  theme_bw(base_size=25)
dev.off()
#dat3 = subset(dat2,distR>1)
