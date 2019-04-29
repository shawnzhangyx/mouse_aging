library(stringr)
statH=commandArgs(trailing=T)[1]
outFile=commandArgs(trailing=T)[2]

#rep=commandArgs(trailing=T)[2]
#if (is.na(rep)) rep = ""
#setwd(paste0("../../analysis/Yang_NMF_method/",tissue,"/",rep))
a=read.table(statH,header=T)
a$barcode = sub(".*\\.(..\\.rep.*)","\\1",a$barcode)
a$stage = substr(a$barcode,1,2)
a$rep = substr(a$barcode,4,7)

a$sample = paste(a$stage,a$rep)

tab = table(a$cluster,a$sample)
tab = sweep(tab,2,colSums(tab),'/')


melted = melt(tab)
g1 = ggplot(melted) + geom_col(aes(factor(Var2),value,fill=factor(Var1)),position="fill")+theme_bw()
g2 = ggplot(melted) + 
  geom_col(aes(factor(Var1),value,fill=factor(substr(Var2,1,3))
    ,color=factor(substr(Var2,4,7))),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"))  +
  theme_bw()

library(gridExtra)
pdf(outFile,height=5,width=10)
g1
g2
#grid.arrange(g1,g2,nrow=1)
dev.off()

