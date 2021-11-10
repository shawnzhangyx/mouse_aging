library(ggplot2)
library(ggrepel)
setwd("/projects/ps-renlab/lamaral/projects/aging_RNA/DH/analysis/repeat_expression/")


lkey = read.table("/projects/ps-renlab/lamaral/software/scTE/mm10.te.bed")
LINES = unique(paste(lkey[which(lkey$V5=="LINE"),4]))
LTR = unique(paste(lkey[which(lkey$V5=="LTR"),4]))
SINE = unique(paste(lkey[which(lkey$V5=="SINE"),4]))
Retroposon = unique(paste(lkey[which(lkey$V5=="Retroposon"),4]))

#key = read.csv("/projects/ps-renlab/lamaral/projects/aging_RNA/DH/analysis/rna_atac_cons_key.csv", header =F)


files = list.files(".", "18.txt")
big = list()
for(f in files) {
  diff_tab = read.table(f)
  cl = gsub("03vs18.txt", "", f)
  cat(f, "\n")
  if ("Ttr" %in% rownames(diff_tab)) {
    diff_tab = diff_tab[-which(rownames(diff_tab)=="Ttr"),]
  }
  diff_tab$celltype = cl
  diff_tab$type = "gene"
  diff_tab[which(rownames(diff_tab) %in% LINES),"type"] = "LINE"   
  diff_tab[which(rownames(diff_tab) %in% Retroposon),"type"] = "Retroposon"   
  diff_tab[which(rownames(diff_tab) %in% SINE),"type"] = "SINE"   
  diff_tab[which(rownames(diff_tab) %in% LTR),"type"] = "LTR"   
  big[[cl]] = diff_tab
}
big_tab = do.call(rbind,big)  

big_tab$type = factor(big_tab$type , levels = c("gene", "LTR", "SINE", "LINE"))
big_tab = big_tab[order(big_tab$type), ]

#big_tab$IAPLTR3 = ""
#big_tab$IAPLTR3[grep("IAPLTR3-int",rownames(big_tab))] = "IAPLTR3-int"
#big_tab$IAPLTR3 = factor(big_tab$IAPLTR3 , levels = c("", "IAPLTR3-int"))

big_tab$gene= sub(".*\\.([^.]*)","\\1",rownames(big_tab))
#big_tab = subset(big_tab, celltype!="Inh.Meis2")
setwd("/mnt/silencer2/home/shz254/projects/mouse_aging/aging_share/figures/repetitive_elements")

big_tab = subset(big_tab,celltype!="NONE")

pdf("DH.repeat_ma_volcano_allclusters.pdf", height = 20, width = 20)

ggplot(big_tab , aes(x=logFC, y=-log10(PValue), color = type)) + geom_point(size = 1, alpha = .5 )  + 
  geom_hline(yintercept=2,linetype="dashed") + 
  geom_label_repel(data=subset(big_tab, PValue<0.01 & type!="gene"), aes(x=logFC,y=-log10(PValue),label=gene)) +
  theme_bw(base_size=25)+theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="plain")) +
  facet_wrap(. ~ celltype, ncol = 5, scales = "free")+ scale_color_manual(values=c("gray", "#E69F00", "#56B4E9", "#f000ff"))

dev.off()

