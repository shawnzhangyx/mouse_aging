infiles = list.files(pattern="bw.tab")

dat = {}
for (f in infiles){
  a = fread(f)
  a$name = f
  dat[[f]] = a
  }

dat2 = data.frame(do.call(rbind,dat))
dat2$pos = as.numeric(sub("s(.*)","\\1",dat2$V1))
dat2$month = sub("(..).metacell_(.*).(..).sorted.rpkm.bw.tab","\\3",dat2$name)
dat2$celltype = sub("(..).metacell_(.*).(..).sorted.rpkm.bw.tab","\\1.\\2",dat2$name)

#agg = aggregate(V6~pos+month,dat2,mean)
 
#ggplot(agg) + geom_line(aes(x=pos,y=V6,group=month,color=month))

avg = aggregate(V6~celltype+month,dat2,mean)

mat = reshape(data=avg,timevar="month",idvar="celltype",direction="wide")

to_vec1 = function(vec){
  vec/sqrt(sum(vec**2))
  }
mat[,2:4] = t(apply(mat[,2:4],1,to_vec1))

dist.mat = dist(mat[,2:4],method="euclidean")

hc = hclust(dist.mat)

mat2 = mat[hc$order,]

melted = melt(mat2)

ggplot(melted) + geom_tile(aes(x=variable,y=factor(celltype,levels=mat2$celltype),
                  ,fill=value))

## normalize a certain 
#mat.norm = reshape(data=dat2,timevar="pos",idvar="celltype",direction="wide")
mat3 = aggregate(V6~celltype,dat2,function(vec){sqrt(sum(vec**2))})
dat2$norm = dat2$V6/mat3$V6[match(dat2$celltype,mat3$celltype)]
dat2$change = dat2$celltype %in% mat2$celltype[1:19]

#read cell type name. 
ctname = read.delim("../celltype_annotation.txt",sep=" ")
ctname$celltype = paste0(ctname$Tissue,".",ctname$Cluster)
ctname$Tname = paste0(ctname$Tissue,".",ctname$Name)

dat2$Tname = ctname$Tname[match(dat2$celltype,ctname$celltype)]
mat2$Tname = ctname$Tname[match(mat2$celltype,ctname$celltype)]


avg.change = aggregate(norm~change+month+pos,dat2,mean)
pdf("normalized_signal_tile.pisd.pdf")
ggplot(dat2) + geom_tile(aes(x=pos,y=factor(Tname,levels=mat2$Tname),
                  ,fill=norm)) + facet_wrap(~month) +
#          scale_fill_gradient2(high="red",mid="white")
                  scale_fill_viridis()


ggplot(avg.change) + geom_line(aes(x=pos,y=norm,group=change,color=change)) + 
    facet_wrap(change~month)

dev.off()
 


