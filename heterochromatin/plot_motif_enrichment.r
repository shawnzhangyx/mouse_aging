#### plot the motif data. 

a =read.delim("DH_FC.UpPeaks.homer.bg/knownResults.txt")

a2 = head(a)
a2$Motif.Name = factor(a2$Motif.Name,levels=rev(a2$Motif.Name))

pdf("enriched_motifs.pdf",height=2,width=8)
ggplot(head(a2)) + geom_col(aes(x=Motif.Name,y=-Log.P.value)) +
  coord_flip()
  dev.off()

