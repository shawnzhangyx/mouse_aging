tissue=commandArgs(trailing=T)[1]
setwd(paste0("../../analysis/snapATAC/",tissue,"/age_diff_edgeR.snap/"))

files = list.files(pattern=".*[up|down].bed",full.names=T)

system("mkdir -p great_chipseq")
output = sub("\\.\\/(.*)\\.bed","./great_chipseq/\\1.great.tsv",files) 
input = gsub("/","%2F",files)
output.scripts = paste(output, " http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=mm10&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2Frenlab.sdsc.edu%2Fyanxiao%2Fmouse_aging%2Fanalysis%2FsnapATAC%2F",tissue,"%2Fage_diff_edgeR.snap%2F",input,sep="")

write.table(output.scripts, "./great_jobs.txt",quote=F,row.names=F,col.names=F)

