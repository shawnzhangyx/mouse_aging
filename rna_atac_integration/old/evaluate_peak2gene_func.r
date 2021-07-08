library(SnapATAC)

load("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/ATAC_imputed.RData")

# correlation analysis
# peaks/DH_peak_83639    chr16   72924125 72925126    #0.6052693

# peak2gene 
# chr16 72810662 72811663 #170.87230

# peak of interest
# chr16 72772036 72773037
# promoter of interest
# peaks/DH_peak_83584 chr16:72307226-72308226

p.x.sp@gmat 

# models <- suppressWarnings(llply(.data=peaks.id, .fun=function(t) summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",], .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE));
#       summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",]
gene.val = p.x.sp@gmat[,which(colnames(p.x.sp@gmat)=="Xkr4")]
data.use = p.x.sp@pmat
t = which(p.x.sp@peak$name =="chr16:72810662-72811663")
t = which(p.x.sp@peak$name =="chr16:72924125-72925126")
t = which(p.x.sp@peak$name =="chr16:72307225-72308226")
t = which(p.x.sp@peak$name =="chr16:72699032-72700033")

models = summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t]>0, x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",]



plotFeatureSingle(p.x.sp, 
        feature.value=p.x.sp@gmat[,which(colnames(p.x.sp@gmat)=="Xkr4")],
        method="umap",
        main="Robo1",
        point.size=0.5, 
        point.shape=19, 
#        down.sample=10000,
        quantiles=c(0.01, 0.99)
 )
## combine RNA and x,y
dat = data.frame(umap1=p.x.sp@umap[,1],umap2=p.x.sp@umap[,2],gene=gene.val,sample=p.x.sp@sample)
dat$age = sub("DH_(..)_rep.","\\1",dat$sample)


ggplot(dat) + geom_point(aes(x=umap1,y=umap2)) + 
  facet_wrap(.~cluster)

dat = data.frame(y=data.use[,t]>0, x=gene.val)
ggplot(dat) + geom_boxplot(aes(x=y,y=x))


comb = data.frame(nRNA = rowSums(p.x.sp@gmat),nATAC=rowSums(p.x.sp@pmat))

