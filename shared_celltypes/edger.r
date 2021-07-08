library(edgeR)

a=read.delim("shared_endo.counts",skip=1)

mat = a[,-c(1:6)]
colnames(mat) = sub("shared.endo.bam.(.*).bam","\\1",colnames(mat))
rownames(mat) = a$Geneid
mat = mat[order(-rowSums(mat)),]

mat2 = mat[-c(1:100),]


y = DGEList(mat2)
y = y[which(rowSums(cpm(y)>1)>=2),]
y = calcNormFactors(y)

d = plotMDS(y)

dat = data.frame(sample=names(d$x),x=d$x,y=d$y)

dat$ct = sub("(..).metacell_(.*).(..).rep.","\\1.\\2",dat$sample)
dat$age = sub("(..).metacell_(.*).(..).rep.","\\3",dat$sample)
dat$rep = sub("(..).metacell_(.*).(..).(rep.)","\\4",dat$sample)

ggplot(dat) + geom_text(aes(x=x,y=y,label=paste0(ct,":",age),color=rep))

