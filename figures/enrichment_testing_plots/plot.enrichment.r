#DG neurons. DH 2 Up
#Muscles. LM 1 Down
#Neutrophils. BM 1

a=read.delim("../../../analysis/snapATAC/DH/age_diff_edgeR.snap/motif.homer.csbg/2.up.homer/knownResults.txt")
a2 = head(a)
a2$Motif.Name = factor(a2$Motif.Name,levels=rev(a2$Motif.Name))

pdf("DH.2.Up.enriched_motifs.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Motif.Name,y=-Log.P.value)) +
  coord_flip() + ggtitle("DG.UP.Motif") +
  theme_bw()
dev.off()

tmp=readLines("../../../analysis/snapATAC/DH/age_diff_edgeR.snap/great_chipseq/2.up.great.red.tsv")
tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)
gobp_red = tab[which(tab$X..Ontology=="GO Biological Process"),]

a2 = gobp_red[which(gobp_red$TotalGenes>5  & gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1 & gobp_red$TotalGenes<1000),]

a2=head(a2)
a2$Desc = factor(a2$Desc, levels=rev(a2$Desc))

pdf("DH.2.Up.enriched_GOBP.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Desc,y=-log10(BinomFdrQ))) +
  coord_flip() + ggtitle("DG.UP.GOBP") +
  theme_bw()
dev.off()


# LM 1 down
a=read.delim("../../../analysis/snapATAC/LM/age_diff_edgeR.snap/motif.homer.csbg/1.down.homer/knownResults.txt")
a2 = head(a)
a2$Motif.Name = factor(a2$Motif.Name,levels=rev(a2$Motif.Name))

pdf("LM.1.Down.enriched_motifs.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Motif.Name,y=-Log.P.value)) +
  coord_flip() + ggtitle("Skm.Down.Motif") +
  theme_bw()
dev.off()

tmp=readLines("../../../analysis/snapATAC/LM/age_diff_edgeR.snap/great_chipseq/1.down.great.red.tsv")
tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)
gobp_red = tab[which(tab$X..Ontology=="GO Biological Process"),]

a2 = gobp_red[which(gobp_red$TotalGenes>5  & gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1 & gobp_red$TotalGenes<1000),]

a2=head(a2)
a2$Desc = factor(a2$Desc, levels=rev(a2$Desc))

pdf("LM.1.Down.enriched_GOBP.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Desc,y=-log10(BinomFdrQ))) +
  coord_flip() + ggtitle("Skm.Down.GOBP") +
  theme_bw()
dev.off()

# BM 1 Up

a=read.delim("../../../analysis/snapATAC/BM/age_diff_edgeR.snap/motif.homer.csbg/1.up.homer/knownResults.txt")
a2 = head(a)
a2$Motif.Name = factor(a2$Motif.Name,levels=rev(a2$Motif.Name))

pdf("BM.1.Up.enriched_motifs.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Motif.Name,y=-Log.P.value)) +
  coord_flip() + ggtitle("Neutrophil.Up.Motif") +
  theme_bw()
dev.off()

tmp=readLines("../../../analysis/snapATAC/BM/age_diff_edgeR.snap/great_chipseq/1.up.great.red.tsv")
tab = read.csv(text=tmp,sep='\t',nrows=length(tmp),stringsAsFactors=F,skip=3)
gobp_red = tab[which(tab$X..Ontology=="GO Biological Process"),]

a2 = gobp_red[which(gobp_red$TotalGenes>5  & gobp_red$RegionFoldEnrich>=1.5 & gobp_red$HyperFdrQ < 0.1 & gobp_red$TotalGenes<1000),]

a2=head(a2)
a2$Desc = factor(a2$Desc, levels=rev(a2$Desc))

pdf("BM.1.Up.enriched_GOBP.pdf",height=3,width=8)
ggplot(head(a2)) + geom_col(aes(x=Desc,y=-log10(BinomFdrQ))) +
  coord_flip() + ggtitle("Neutrophil.Up.GOBP") +
  theme_bw()
dev.off()

