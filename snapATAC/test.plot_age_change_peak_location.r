setwd("../../analysis/snapATAC/DH")

files = list.files(path="age_diff_edgeR.snap",pattern="edger.txt",full.names=T)

file = files[6]

chr.len = read.table("/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len = rbind(chr.len,chr.len)
chr.len$logFC = rep(c(-1,1),each=nrow(chr.len)/2)
colnames(chr.len)[1] = "chr"


diff_dict = list()

pdf("diff_peaks.chr_location.pdf")

for (file in files){
print(file)

a=read.csv(file) #"age_diff_edgeR.snap/2.down.bed")
a$chr = sub("(.*):(.*)-(.*)","\\1",a$X)
a = a[which(a$PValue<0.001 & a$chr %in% chr.len$chr),]

if ( nrow(a) > 0) {
a$sample= file 
diff_dict[[file]]  = a 

a$pos = floor(as.integer(sub("(.*):(.*)-(.*)","\\2",a$X))/1e6)

b = a%>% count(chr,pos,logFC>0)
colnames(b)[3] = "logFC"

print(
ggplot(b) + 
  geom_rect(data=chr.len, aes(xmin=0,ymin=0,xmax=V2/1e6,ymax=5),fill=NA,color="black") + 
  geom_col(data=b,aes(x=pos,y=n)) + 
  facet_grid(chr~logFC>0) + 
  theme(
    panel.background = element_rect(fill = NA, colour = NA),
    panel.grid = element_blank(),
#    axis.text.x = element_blank(),
    strip.text.y = element_text(angle = 0),
#    axis.ticks.x=element_blank()
  ) + 
  ggtitle(file) 
)

}

}

dev.off()


out = do.call(rbind,diff_dict)
out$pos = floor(as.integer(sub("(.*):(.*)-(.*)","\\2",out$X))/1e6)

b = out %>% count(chr,pos,logFC>0)

ggplot(out) +geom_bar(aes(x=pos,fill=sample),stat="count",position="stack") +
 facet_grid(chr~logFC>0) 



