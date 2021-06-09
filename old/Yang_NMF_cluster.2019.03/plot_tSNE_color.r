tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/Yang_NMF_method/",tissue))
a=read.table(paste0(tissue,".statH"))

a$V3 = a$V3+1
a$stage = substr(a$V1,1,2)
a$rep = substr(a$V1,4,7)

a$sample = paste(a$stage,a$rep)


#b = read.table(paste0(tissue,".tsne.xy"))
b = read.table(paste0(tissue,".umap.xy"))
colnames(b) = c("x","y")
a=cbind(a,b)

png("tSNE_plot_color_by_age.png",height=640,width=1280)
ggplot(a) + geom_point(aes(x,y,color=stage),size=1) +
  facet_wrap(.~stage) +
  theme_bw()
dev.off()
png("tSNE_plot_color_by_rep.png",height=640,width=1280)
ggplot(a) + geom_point(aes(x,y,color=rep),size=1) +
  facet_wrap(.~rep) +
  theme_bw() 
dev.off()

