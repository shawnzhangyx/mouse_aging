setwd("../../analysis/snapATAC/DH")

files = list.files(path="PePr",pattern="bed",full.names=T)

file = files[6]

chr.len = read.table("/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len = rbind(chr.len,chr.len)
chr.len$logFC = rep(c(-1,1),each=nrow(chr.len)/2)
colnames(chr.len)[1] = "chr"


diff_dict = list()

pdf("diff_peaks.chr_location.pdf")

for (file in files){
print(file)

a=read.table(file) #"age_diff_edgeR.snap/2.down.bed")
a$chr = a$V1
#a = a[which(a$PValue<0.001 & a$chr %in% chr.len$chr),]

a$sample= file 
diff_dict[[file]]  = a 

a$pos = floor(a$V2/1e6)

b = a%>% count(chr,pos)
#colnames(b)[3] = "logFC"

print(
ggplot(b) + 
  geom_rect(data=chr.len, aes(xmin=0,ymin=0,xmax=V2/1e6,ymax=5),fill=NA,color="black") + 
  geom_col(data=b,aes(x=pos,y=n)) + 
  facet_grid(chr~.) + 
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


dev.off()


out = do.call(rbind,diff_dict)
out$chr = as.character(out$V1)
out$pos = floor(out$V2/1e6)
out$chip1 = grepl("chip1",out$sample)
out$cluster = sub(".*metacell_(.*).03vs18_.*","\\1",out$sample)

out = out[which(out$chr %in% paste0("chr",c(1:19,"X"))),]
out = out[which(out$cluster !="10"),]
b = out %>% count(chr,pos,chip1)

ggplot(out) +geom_bar(aes(x=pos,fill=cluster),stat="count",position="stack") +
 facet_grid(chr~chip1) +
 theme(
     panel.background = element_rect(fill = NA, colour = NA),
     panel.grid = element_blank(),
     strip.text.y = element_text(angle = 0)
             )



