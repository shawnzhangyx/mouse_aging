tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/Yang_NMF_method/",tissue))


a=read.delim("heart.statH")
p1 = read.delim("rep1/heart.statH")
p2 = read.delim("rep2/heart.statH")
a$class0 = a$class0+1
p1$class0 = p1$class0+1
p2$class0 = p2$class0+1

p1$pooled.rank = a$class0[match(p1$X..xgi,a$X..xgi)]
table(p1$class0,p1$pooled.rank)
#p1$corrected_class = c(2,1,3,4,5)[p2$class0]
p1$corrected_class  = p1$class0
table(p1$corrected_class,p1$pooled.rank)


p2$pooled.rank = a$class0[match(p2$X..xgi,a$X..xgi)]
table(p2$class0,p2$pooled.rank)
p1$corrected_class = c(2,1,3,4,5)[p1$class0]
table(p1$corrected_class,p1$pooled.rank)


table(a$class0)/nrow(a)
table(p1$class0)/nrow(p1)
table(p2$class0)/nrow(p2)

