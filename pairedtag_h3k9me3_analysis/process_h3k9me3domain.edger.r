library(edgeR)

a=read.delim("SICER-W5000-G10000.reps.read.counts",skip=1)
samples = sub("bam_merge.(.*).rep..bam","\\1",colnames(a)[-c(1:6)])
cts = sub("(.*)...","\\1",samples)

# remove chrY
a = a[which(a$Chr!= "chrY"),]
#a = a[which(a$Length>100000),]
ct = "ASC"

for (ct in c("ASC","CA1","CA23","Claustrum","DG","InNeuPvalb","InNeuSst","InNeuVip","L23","L4","L5Deptor","L5Fezf2","L5Parm1","L5Unknown","L6","OGC")){
print(ct)
idx = which(cts ==ct)
cnts2 = a[,idx+6]
row.names(cnts2) = paste0(a$Chr,":",a$Start,"-",a$End)

y = DGEList(cnts2)
#y$sample$lib.size=total[idx]
#y = y[which(rowSums(cpm(y)>1)>=2),]
#y = calcNormFactors(y)


groups = samples[idx]
design = model.matrix(~0+groups)
y<-estimateCommonDisp(y)
y<-estimateGLMTagwiseDisp(y,design)
fit_tag = glmFit(y,design)

lrt = glmLRT(fit_tag, contrast =c(-1,1))
fdr = p.adjust(lrt$table$PValue,method="BH")
cpms = cpm(y) 
cpms = (cpms[,c(1,3)] + cpms[,c(2,4)])/2
colnames(cpms) = c("M03","M18")
out = cbind(cpms,lrt$table,fdr)
out = out[which(a$Length > 100000),]
out = out[order(out$PValue),]

write.csv(out,paste0("edger/",ct,".peaks.age_diff.csv"))
}

down=0
up = 0
out2 = list()
for (ct in c("ASC","CA1","CA23","Claustrum","DG","InNeuPvalb","InNeuSst","InNeuVip","L23","L4","L5Deptor","L5Fezf2","L5Parm1","L5Unknown","L6","OGC")){
#print(ct)
out = read.csv(paste0("edger/",ct,".peaks.age_diff.csv"))
print(c(ct,length(which(out$logFC< -0.6 & out$fdr<0.2)), length(which(out$logFC>0.6 & out$fdr<0.2))))
#down=down+length(which(out$logFC< -0.6 & out$fdr<0.2))
#up = up + length(which(out$logFC>0.6 & out$fdr<0.2))
out$ct = ct 
out$ct2 = ifelse(ct %in% c("ASC","OGC"), "Glia",ifelse(ct %in% c("InNeuPvalb","InNeuSst","InNeuVip"), "InN","ExN"))

out2[[ct]] = out
}
dat = do.call(rbind,out2)
dat$start = as.numeric(sub("chr.*:(.*)-(.*)","\\1",dat$X))
dat$end = as.numeric(sub("chr.*:(.*)-(.*)","\\2",dat$X))
dat$size = log10(dat$end-dat$start)

pdf("h3k9me3_domain.scatter.pdf",height=6,width=8)
ggplot(dat) + geom_point(aes(x=logFC,y=-log10(fdr),size=size),color="grey",alpha=0.6) + 
#  geom_point(data=subset(dat,fdr<0.05 ),aes(x=logFC,y=-log10(fdr),color=ct2)) + scale_y_continuous(limits=c(0,20)) +
 geom_point(data=subset(dat,fdr<0.05 & ( logFC>0.6 | logFC < -0.6)),aes(x=logFC,y=-log10(fdr),color=ct2,size=size),alpha=0.6) + scale_y_continuous(limits=c(0,20)) + 
  xlim(-4,4) + geom_vline(xintercept=0.1,linetype="dashed") + 
  theme_classic(base_size=28)
  dev.off()



sig = dat[which(dat$fdr<0.05 & abs(dat$logFC) > 0.6),]

sig.uniq = sig[!duplicated(sig[,c("X","ct2")]),]






