setwd("../../analysis/snapATAC/FC/snapFiles/")

library(SnapATAC)
load("FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.cluster.RData")
load("FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData")

m1= x.sp@smat@dmat

x.after.sp@cluster

tab  = data.frame(x.sp@sample,x.sp@barcode,x.after.sp@cluster,x.sp@umap, x.after.sp@umap)
tab$age = substr(tab$x.sp.sample,4,5)
tab$rep = substr(tab$x.sp.sample,7,10)

ggplot(subset(tab,x.after.sp.cluster==1)) + 
  geom_point(aes(x=umap.1,y=umap.2,color=age),alpha=0.5) +
  xlim(0,8) + ylim(-10,-5)

ggplot(subset(tab,x.after.sp.cluster==1)) + 
  geom_density_2d_filled(aes(x=umap.1,y=umap.2)) +
  facet_grid(age~.) +
#  scale_fill_manual(values=hsv(1, seq(0,1,length.out = 11) , 1)) +
  #geom_point(aes(x=umap.1,y=umap.2),alpha=0.1) + 
  xlim(0,8) + ylim(-10,-5)

ggplot(subset(tab,x.after.sp.cluster==1)) +
#  geom_point(aes(x=umap.1,y=umap.2),alpha=0.1) +
  geom_density_2d_filled(aes(x=umap.1,y=umap.2)) +
  facet_grid(age~rep) +
  xlim(0,8) + ylim(-10,-5)




ggplot(subset(tab,x.after.sp.cluster==1)) +
  geom_point(aes(x=umap.1.1,y=umap.2.1,color=age),alpha=0.5) +
  xlim(-3,3) + ylim(2,8)

ggplot(subset(tab,x.after.sp.cluster==1)) +
   geom_point(aes(x=umap.1.1,y=umap.2.1),alpha=0.1) +
  geom_density_2d_filled(aes(x=umap.1.1,y=umap.2.1)) +
  #scale_fill_manual(values=hsv(1, seq(0,1,length.out = 11) , 1)) +
  facet_grid(age~rep) +
  xlim(-3,3) + ylim(2,8)


samples = x.sp@sample[which(x.after.sp@cluster==1)]


m1= x.sp@smat@dmat[which(x.after.sp@cluster==1),]

m2= x.after.sp@smat@dmat[which(x.after.sp@cluster==1),]

ct_list = list()
for (sample in unique(samples)) {
  ct_list[[sample]] = apply(m1[which(samples==sample),],2,median)

}

out = do.call(rbind,ct_list)


dist(out) 






