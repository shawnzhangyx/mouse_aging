tissue =  commandArgs(trailing=T)[1]
rank = commandArgs(trailing=T)[2]

genesets =  c("GO Biological Process","Mouse Phenotype","MGI Expression: Detected")

setwd(paste0("../../analysis/Yang_NMF_method/",tissue,"/R",rank,"/age_diff_bycluster/great"))

fname = "dorsal_hippocampus.C3.peaks.up.great.red.tsv"

dat.list = {}
for (fname in list.files(pattern="red.tsv")) {
file = readLines(fname)
file = file[-c(1:4,(length(file)-21):length(file))]
tab = read.csv(text=file,sep='\t',nrows=length(file),stringsAsFactors=F)
tab = tab[which(tab$X..Ontology %in% genesets),]


gobp = tab[which(tab$X..Ontology=="GO Biological Process"),]
pheno = tab[which(tab$X..Ontology=="Mouse Phenotype"),]

colnames(gobp)[c(2,20)] = colnames(pheno)[c(2,20)] = c("Geneset.ID","N.genes")
gobp_red = gobp
gobp_red = gobp_red[which(gobp_red$N.genes>5  & gobp_red$RegionFoldEnrich>=2),]
pheno = pheno[which(pheno$N.genes>5 &pheno$RegionFoldEnrich >=2),]   # & pheno$HyperFdrQ < 0.05),]

tab = data.frame(fname,gobp_red[1:5,c("Desc","BinomP")],pheno[1:5,c("Desc","BinomP")])
dat.list[[fname]] = tab

}

dat = do.call(rbind,dat.list)

dat$fname = sub(".*C(.*).peaks.(.*).great.red.tsv","\\1.\\2",dat$fname)

dat2 = dat[,1:3]
dat2$BinomP = -log10(dat2$BinomP)

library(grid)
library(gridExtra)
pdf("great.table.pdf",height=20,width=20)
grid.table(dat2)
dev.off()
