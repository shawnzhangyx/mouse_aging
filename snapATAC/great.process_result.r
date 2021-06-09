tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/"))

files = list.files(path="great_chipseq",pattern = "*great.red.tsv",full.names=T,recursive=T)

dat = list()

for (file in files){
  print(file)
  cluster = sub("great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\1",file)
  category = sub("great_chipseq\\/(.*)\\.(.*).great.red.tsv","\\2",file)
  tmp=readLines(file)
  tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)

  gobp = tab[which(tab$X..Ontology=="GO Biological Process"),]
  if (length(gobp) <1) { next }
  gobp_red = gobp
  gobp_red = gobp_red[which(gobp_red$TotalGenes>5  & gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1),]

  tmp = gobp_red[1:5,c("ID","Desc","BinomFdrQ")]

  dat[[file]] = data.frame(cluster=cluster,category=category,tmp)

  }

out = do.call(rbind,dat)
write.csv(out,"DE_Peaks.great_results.csv")

