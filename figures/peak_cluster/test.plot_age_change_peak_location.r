#setwd("../../analysis/snapATAC/DH")

files = list.files(path="/projects/ps-renlab/yanxiao/projects/mouse_aging/aging_share/clustered_diff_peak_locations/",pattern="smth20",recursive=T,full.names=T)
files = files[!grepl("old_May5",files)]
files = files[!grepl("old_2020May6",files)]
files = files[grepl("(DH|FC)",files)]


file = files[2]
chr.len = read.table("/projects/ps-renlab/share/bwa_indices/mm10.fa.fai")
chr.len$color = c(rep(c("black","grey"),10),"black")
colnames(chr.len)[1] = "chr"


chr.window = NULL
size = 10^5
for (i in 1:nrow(chr.len)){
    chr.window = c(chr.window, paste0(chr.len[i,1],":",0:ceiling(chr.len[i,2]/size) ))
    }

diff_dict = list()


for (file in files){
print(file)
a=data.frame(fread(file))
a$pos = factor(paste0(a$chr,":",a$window),levels=chr.window)
#a$tissue = sub(".*\\/(..)\\/smth20_scores_table_(.*).txt","\\1",file)
#a$direct = sub(".*\\/(..)\\/smth20_scores_table_(.*).txt","\\2",file)
#a$color = chr.len$color[match(a$chr,chr.len$chr)]
a$color = ifelse(grepl("Up",file),"red","blue")
name = sub(".*\\/(..)\\/smth20_scores_table_(.*).txt","\\1.\\2",file)

diff_dict[[name]] = a
}

out = do.call(rbind,diff_dict)
out = out[!is.na(out$pos),]
out$posN = as.numeric(out$pos)
out$chr = factor(out$chr, levels=paste0("chr",c(1:19,"X","Y")))


sig = read.delim("clustered_diff_peak_all.smth20_gt0.4.size_gt100k.txt")
sig$chr = sub("(.*):(.*)-(.*)","\\1",sig$locations)
sig$pos = factor(paste0(sig$chr,":",ceiling( (sig$start+sig$end)/2/size)),levels=chr.window)
sig$posN = as.numeric(sig$pos)
sig$chr = factor(sig$chr, levels=paste0("chr",c(1:19,"X","Y")))

## plot H3K9me3 domain. 
k9 = read.table("/mnt/silencer2/home/shz254/projects/mouse_aging/data/encode/histone/peaks/FB.P0.H3K9me3-W1000-G3000-FDR0.01-island.bed")
colnames(k9) = c("chr","start","end","score")
k9$start_pos = factor(paste0(k9$chr,":",ceiling(k9$start/size)),levels=chr.window)
k9$end_pos = factor(paste0(k9$chr,":",ceiling(k9$end/size)),levels=chr.window)
k9$size = k9$end-k9$start
k9$start_posN = as.numeric(k9$start_pos)
k9$end_posN = as.numeric(k9$end_pos)
k9 = k9[which(k9$size>=100e3),]
k9$chr = factor(k9$chr, levels=paste0("chr",c(1:19,"X","Y")))
# write the large k9 domains. 
write.table(k9[,c("chr","start","end","score","size")],"FB.P0.H3K9me3.gt100k_domain.bed",row.names=F,col.names=F,quote=F,sep="\t")


#png("DH_FC.diff_peaks.chr_location.png",units="in",width=16,height=6,res=300)
pdf("DH_FC.diff_peaks.chr_location.pdf",height=6,width=16)
ggplot(out) +
  geom_line(aes(x=posN,y=smth20_score,color=color),alpha=0.5) + 
  scale_color_manual(values=c("blue","red")) + 
  geom_jitter(data=subset(sig,direct=="Down"),aes(x=posN,y=-0.25),shape=25,color="blue",fill="blue",width=0,height=0.05) +
  geom_jitter(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.25),shape=24,color="red",fill="red",width=0,height=0.05) +
#  geom_text(data=subset(sig,direct=="Up"),aes(x=posN,y=-0.35,label=pos),angle=45,color="red") +
  geom_segment(data=k9,aes(x=start_posN,y=-0.15,xend=end_posN,yend=-0.15),size=5) +
  geom_hline(yintercept=0.4,linetype="dashed") +
  ylim(-0.5,2) +
  theme_bw() + 
  theme( axis.ticks.x= element_blank(),
  axis.text.x = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.grid = element_blank(),
  panel.border = element_rect(colour="grey")
  )+
  facet_grid(.~chr,scales="free_x",space="free_x") 
dev.off()

#out = out[which(out$chr=="chr9"),]
#sig = sig[which(sig$chr=="chr9"),]
#k9 = k9[which(k9$chr=="chr9"),]
#k9$start_pos = factor(paste0(k9$chr,":",ceiling(k9$start/size)),levels=chr.window)
#k9$end_pos = factor(paste0(k9$chr,":",ceiling(k9$end/size)),levels=chr.window)


#ggplot(out) +geom_line(aes(x=pos,y=smth20_score,color=color)) +
#  geom_point(data=subset(sig,direct=="Down"),aes(x=pos,y=-0.25),shape=25,color="blue",fill="blue") +
#  geom_point(data=subset(sig,direct=="Up"),aes(x=pos,y=-0.25),shape=24,color="red",fill="red") +
#  geom_segment(data=k9,aes(x=start_pos,y=-0.5,xend=end_pos,yend=-0.5),size=5)





#ggplot(sig) +geom_point(aes(x=as.numeric(pos),y=-1),shape=2)
#ggplot(sig) +geom_point(aes(x=as.numeric(pos),y=paste0(tissue,clust),color=tissue))


