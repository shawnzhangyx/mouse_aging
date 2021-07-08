a=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/pearson_corr_concat.txt"))

b=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/spearman_corr_concat.txt"))

d=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/weighted_regression.txt"))


a = a[order(a$V10),]
b = b[order(b$V10),]
d = d[order(d$Pvalue),]

a$fdr = p.adjust(a$V10,method="BH")
b$fdr = p.adjust(b$V10,method="BH")
d$fdr = p.adjust(d$Pvalue,method="BH")

table(a$fdr<0.05)
table(b$fdr<0.05)
table(d$fdr<0.05)


robo1.a = a[which(a$V1=="Robo1"),]
robo1.b = b[which(b$V1=="Robo1"),]
robo1.d = d[which(d$Geneid=="Robo1"),]

robo1.a[grep("peaks/DH_peak_83584",robo1.a$V5),]
robo1.b[grep("peaks/DH_peak_83584",robo1.b$V5),]
robo1.d[grep("peaks/DH_peak_83584",robo1.d$Peakid),]


robo1.d[which(robo1.d$ATAC_start==72772036),]
robo1.a[which(robo1.a$V7==72772036),]
robo1.b[which(robo1.b$V7==72772036),]

robo1.d$Start =NULL
robo1.d$End =NULL

#72772036
