library(amap)
library(edgeR)
set.seed(1)

mark="all_celltypes.thresh.filter"
setwd(paste0("../../analysis/all_celltypes.all_tissue"))
a=data.frame(fread(paste0(mark,".counts"),skip=1))
#b=data.frame(fread(paste0(mark,"_ATACpeak.overlaped.bed")))
cnts = a[,-c(1:6)]
grps = paste0("C",sub(paste0("X(.*).metacell.bam"),"\\1",colnames(cnts)))
colnames(cnts) = grps
rownames(cnts) = paste0(a$Chr,":",a$Start,"-",a$End)

y <- DGEList(counts=cnts);
y = calcNormFactors(y)

cpms = cpm(y)
cpms = cpms[,order(as.numeric(sub("C(.*)","\\1",colnames(cpms))))]
cpms = cpms[,-29] # remove the bad clustere. 

norm = log10(cpms+1e-2)
norm[norm<0] = 0
norm[norm>2] = 2
norm = norm[which(rowSums(norm)>0),]

#library(doParallel)
#registerDoParallel(cores=50)
#withinss = foreach(i=2:50) %dopar% {
#print(i)
#km = kmeans(x=norm,centers=i,iter.max=100)
#km$tot.withinss
#}

#totss = kmeans(x=norm,centers=1,iter.max=100)$totss
#dat = data.frame(k = 2:50,wss=do.call(rbind,withinss),tss=totss)
#pdf(paste0(mark,".raw.ss.pdf"))
#ggplot(dat) + geom_point(aes(x=k,y=wss/tss))
#dev.off()
#write.csv(dat,paste0(mark,".cluster.ss.csv"))

## pick 35 clusters.
km = kmeans(x=norm,centers=35,iter.max=100)

write.csv(cbind(norm,km$cluster),paste0(mark,".kmeans.cluster.csv"))


#km = Kmeans(x=norm,centers=30,iter.max=20,method="correlation")
hc = hclust(as.dist(1-cor(norm)))

norm2_sorted = norm[order(km$cluster),hc$order]
cluster_sorted = sort(km$cluster)
spike_list = NULL
## find spike method 2
for (idx in 1:length(unique(km$cluster))){
 temp = norm2_sorted[which(cluster_sorted==idx),]
 tab = table(apply(temp,1,which.max))
 spike = as.numeric(names(tab)[which.max(tab)])
 spike_list = c(spike_list,spike)
}


cluster_ordered = rank(spike_list)[match(km$cluster,1:length(unique(km$cluster)))]

norm_ordered = norm[order(cluster_ordered),]
norm_ordered2 = norm_ordered[seq(1,nrow(norm_ordered),10),]

melted = melt(norm_ordered)
melted$cl = km$cluster[order(cluster_ordered)]
## reorganize the clusters manually. 
melted$cl = factor(melted$cl, levels=c(26,4,34,5,10,35,9,3,22,25,29,13,24,27,32,21,31,23,30,2,17,33,1,14,16,8,11,12,19,28,18,20,7,15,6))
melted = melted[order(melted$cl),]

melted$Var1 = factor(melted$Var1,levels=unique(melted$Var1))
melted$Var2 = factor(melted$Var2,levels=paste0("C",c(1:28,30:33))[hc$order] )
#melted$cl = factor(melted$cl, levels=c(1:35)[order(spike_list)])

#system("mkdir H3K9me3_clusters")
library(viridis)
png(paste0(mark,"_ATACpeak_clusters.log10.png"),units="in",height=25,width=30,res=300)
#pdf(paste0(mark,"_ATACpeak_clusters.log10.pdf"),height=10,width=12)
ggplot(melted) +geom_tile(aes(x=Var1,y=Var2,fill=value))+
  scale_fill_viridis() +
  geom_point(aes(x=Var1,y=-1,color=factor(cl)))+
  facet_grid(.~cl,scales="free_x",space="free_x") +
  theme_bw(base_size=20) + 
  theme(axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.spacing = unit(0, "lines"))
dev.off()


avg = aggregate(value~Var2+cl,melted,median)
avg$value[avg$value<0] = 0

pdf("avg.log10.pdf")
ggplot(avg) +geom_tile(aes(x=cl,y=Var2,fill=value))+
  scale_fill_viridis() +
  geom_point(aes(x=cl,y=-1,color=factor(cl)))#+
#  theme(axis.text.x = element_blank(),
#  axis.ticks.x = element_blank())
dev.off()

