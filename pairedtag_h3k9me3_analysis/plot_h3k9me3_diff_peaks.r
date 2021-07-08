
a1=read.csv("edger/L23.peaks.age_diff.csv")
a1$ct = "L23"
a2=read.csv("edger/L4.peaks.age_diff.csv")
a2$ct = "L4"
a3=read.csv("edger/ASC.peaks.age_diff.csv")
a3$ct = "ASC"
a = rbind(a1,a2,a3)

a$chr = sub("(.*):(.*)-(.*)","\\1",a$X)
a$start = as.numeric(sub("(.*):(.*)-(.*)","\\2",a$X))
a$end = as.numeric(sub("(.*):(.*)-(.*)","\\3",a$X))

#a2 = 

ggplot(subset(a,fdr<0.05)) + geom_rect(aes(xmin=start,ymin=0,xmax=end,ymax=logFC,fill=logFC),size=5) +
  scale_fill_gradient2() + 
  facet_grid(ct~chr) 


