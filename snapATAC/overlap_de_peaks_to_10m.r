setwd("../../analysis/snapATAC/FC")

cl = 1

ov_list = NULL
random_list = NULL
for(cl in 1:17){
a =data.frame(fread(paste0("age_diff_edgeR.snap/",cl,".edger.txt")))

b1 =data.frame(fread(paste0("age_diff_edgeR.snap.03m.vs.10m/",cl,".edger.txt")))
b2 =data.frame(fread(paste0("age_diff_edgeR.snap.10m.vs.18m/",cl,".edger.txt")))

t = length(which(a$PValue<0.01))
#random sample t peaks. 
r = sample(a$V1,t)
#o1 = length(which(a$V1[which(a$PValue<0.01)] %in% b1$V1[which(b1$PValue<0.05)]))
#o2 = length(which(a$V1[which(a$PValue<0.01)] %in% b2$V1[which(b2$PValue<0.05)]))
o =  length(which( a$V1[which(a$PValue<0.01)] %in% c(b1$V1[which(b1$PValue<0.05)],b2$V1[which(b2$PValue<0.05)])))
ct = length(which(r %in% c(b1$V1[which(b1$PValue<0.05)],b2$V1[which(b2$PValue<0.05)])))
#print(paste(cl, (o1+o2)/t ))
print(paste(cl, o/t, ct/t ))

ov_list = c(ov_list,o/t)
random_list = c(random_list,ct/t)
}

cl = 1

a =data.frame(fread(paste0("age_diff_edgeR.snap/",cl,".edger.txt")))
b1 =data.frame(fread(paste0("age_diff_edgeR.snap.03m.vs.10m/",cl,".edger.txt")))
b2 =data.frame(fread(paste0("age_diff_edgeR.snap.10m.vs.18m/",cl,".edger.txt")))


t0 = a[which(a$PValue<0.01),]
t1 = b1[which(b1$V1 %in% a$V1[which(a$PValue<0.01)]),]

t2 = b2[which(b2$V1 %in% a$V1[which(a$PValue<0.01)]),]

m = merge(t0[,c(1,8)],t1[,c(1,8)],by="V1")
m = merge(m, t2[,c(1,8)],by="V1")


ggplot(m) + geom_point(aes(x=logFC.x,y=logFC.y))
ggplot(m) + geom_point(aes(x=logFC.x,y=logFC))
ggplot(m) + geom_point(aes(x=logFC.y,y=logFC,color=logFC.x>0))

