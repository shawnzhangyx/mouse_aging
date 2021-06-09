setwd("../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/")
summ2 = read.delim("age_diff_bycluster/dorsal_hippocampus.peaks.age_diff.summary.txt")

#labels_o = rev(c("C10","C15","C6","C2","C5","C14","C12","C3","C7","C13",
#             "C11","C4","C1","C8","C9"))
labels_o = rev(c("C10","C15","C6","C2","C14","C5","C12","C3","C7","C9",
                   "C1","C8","C11","C13","C4"))

summ2$cluster = factor(summ2$cluster,levels=labels_o)
summ2 = summ2[order(summ2$cluster),]

summ3 = gather(summ2,Change,Num,Down_Reg,Up_Reg,)

pdf("age_diff_bycluster/dorsal_hippocampus.peak.age_diff.summary.pdf",width=3,height=4)
ggplot(summ3) + geom_col(aes(x=cluster,y=Num,fill=Change),position="stack") + 
    coord_flip() + theme_bw() #+ scale_fill_brewer(palette="Spectral")
dev.off()


setwd("../../analysis/Yang_NMF_method/frontal_cortex/R15/")
summ2 = read.delim("age_diff_bycluster/frontal_cortex.peaks.age_diff.summary.txt")


labels_o = rev(c("C12","C4","C2","C5","C6","C3","C7","C9","C11","C8",
                   "C13","C1","C10","C14","C15"))
summ2$cluster = factor(summ2$cluster,levels=labels_o)
summ2 = summ2[order(summ2$cluster),]

summ3 = gather(summ2,Change,Num,Down_Reg,Up_Reg,)

pdf("age_diff_bycluster/frontal_cortex.peak.age_diff.summary.pdf",width=3,height=4)
ggplot(summ3) + geom_col(aes(x=cluster,y=Num,fill=Change),position="stack") +
    coord_flip() + theme_bw() #+ scale_fill_brewer(palette="Spectral")
    dev.off()

