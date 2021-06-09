library(ggrepel)
dir=commandArgs(trailing=T)[1]
#output=commandArgs(trailing=T)[2]

setwd(dir)
files=list.files(pattern=".txt")
#file = files[1]
for (file in files) {
print(file)  
a =read_delim(file,col_names=F,delim=" ")
b = a %>% group_by(X3) %>% summarize(median(X1),median(X2)) %>% rename("X1"="median(X1)","X2"="median(X2)")

svg(sub(".txt",".svg",file))
print(ggplot(a) + geom_point(aes(x=X1,y=X2,color=factor(X3+1)),size=0.5) +
    geom_label_repel(data=b,aes(x=X1,y=X2,label=X3+1,color=factor(X3+1))) +
    scale_color_discrete(name="cluster") +
    ggtitle(file) +
    theme( legend.position="none")
)          
dev.off()

png(sub(".txt",".png",file))
print(ggplot(a) + geom_point(aes(x=X1,y=X2,color=factor(X3+1)),size=0.5) +
#    geom_label_repel(data=b,aes(x=X1,y=X2,label=X3+1,color=factor(X3+1))) +
    scale_color_discrete(name="cluster") +
    ggtitle(file) + xlab("UMAP-1") + ylab("UMAP-2") + 
    theme_void() +
    theme( legend.position="none")
)
dev.off()

}


