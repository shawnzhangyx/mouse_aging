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
  if (nrow(tmp)>0){
  dat[[file]] = data.frame(cluster=cluster,category=category,tmp)
  }
  }

out = do.call(rbind,dat)
write.csv(out,"DE_Peaks.great_results.csv")


# plot it. 
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
  tmp = gobp_red[,c("ID","Desc","BinomFdrQ")]

  if (nrow(tmp)>0){
    tmp$rank = 1:nrow(tmp)
  dat[[file]] = data.frame(cluster=cluster,category=category,tmp)
  }
  }

out = do.call(rbind,dat)

library(ggrepel)
pdf(paste0(tissue,".great_results.pdf"),height=50,width=10)
ggplot(out)  + geom_line(aes(x=rank,y=-log10(BinomFdrQ))) + 
  geom_point(aes(x=rank,y=-log10(BinomFdrQ))) +
  geom_label_repel(data=subset(out,rank <=5),aes(x=rank,y=-log10(BinomFdrQ),label=Desc),
 #   nudge_x = .15,
   box.padding = 0.5,
   force = 10,
#    nudge_y = 1,
#    max.overlaps= Inf,
    ) + 
  facet_grid(cluster~category,scales="free_y")
dev.off()


