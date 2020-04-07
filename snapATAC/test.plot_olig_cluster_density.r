
# x -5 0
# y 2 7 
BIN = 100
mat1 = mat2=  matrix(0,nrow=BIN,ncol=BIN)
xseq = seq(-5,0,5/BIN)
yseq = seq(2,7,5/BIN)


tmp = b[which(b$stage==3),]
for (x in 1:(BIN)){
  print(x)
  for (y in 1:(BIN)){
    mat1[x,y] = length(which(tmp$umap.1>=xseq[x] & tmp$umap.1<xseq[x+1] &
      tmp$umap.2>=yseq[y] & tmp$umap.2< yseq[y+1]))
  }
  }


tmp = b[which(b$stage==18),]

for (x in 1:(BIN)){
  print(x)
  for (y in 1:(BIN)){
    mat2[x,y] = length(which(tmp$umap.1>=xseq[x] & tmp$umap.1<xseq[x+1] &
      tmp$umap.2>=yseq[y] & tmp$umap.2< yseq[y+1]))
  }
  }


dat1 = melt(mat1)
dat2 = melt(mat2)

dat = merge(dat1,dat2,by=c("Var1","Var2"))


dat$d1 = dat$value.x/sum(dat$value.x)
dat$d2 = dat$value.y/sum(dat$value.y)


g1 = ggplot(dat)+ geom_tile(aes(x=Var1,y=Var2,fill=d1)) +
  scale_fill_gradient(low='white',high='red')
g2 = ggplot(dat)+ geom_tile(aes(x=Var1,y=Var2,fill=d2)) + 
  scale_fill_gradient(low='white',high='red')

#library(gridExtra)
grid.arrange(g1,g2)


ggplot(dat)+ geom_tile(aes(x=Var1,y=Var2,fill=d1-d2)) +
  scale_fill_gradient2(low='blue',mid="grey",high='red')

