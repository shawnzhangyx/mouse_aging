setwd('../../analysis/repeats_RepEnrich2/')

a = read.table("frontal_cortex_3m/frontal_cortex_3m_fraction_counts.txt")
a2 =read.table("frontal_cortex_18m/frontal_cortex_18m_fraction_counts.txt")

a$V5 = a2$V4
a = a[-which(a$V2 == "Simple_repeat"),]
a$V4 = a$V4/sum(a$V4)* 1e6
a$V5 = a$V5/sum(a$V5)* 1e6
a$log2FC = log2((a$V5+1)/(a$V4+1)) 


library(plotly) 

library(ggrepel)
p = ggplot(a) + geom_point( aes(x=V4,y=V5,color=V3,label=V1)) + 
  geom_abline(intercept=0,slope=1,color='grey') +
  scale_x_log10() +scale_y_log10()

ggplotly(p)

