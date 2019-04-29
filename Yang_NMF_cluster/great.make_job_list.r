path = commandArgs(trailing=T)[1] 
tissue =  commandArgs(trailing=T)[2]
rank = commandArgs(trailing=T)[3]

files = list.files(path=path,pattern="(up|down).bed",recursive=T)
output = paste0("great/",sub("bed","great.tsv",files))
input = gsub("/","%2F",files)
output.scripts = paste0(path,"/",output, " http://bejerano.stanford.edu/great/public/cgi-bin/greatStart.php?outputType=batch&requestSpecies=mm10&requestName=Example+Data&requestSender=Client+A&requestURL=http%3A%2F%2Frenlab.sdsc.edu%2Fyanxiao%2Fmouse_aging%2Fanalysis%2FYang_NMF_method%2F",tissue,"%2FR",rank,"%2Fage_diff_bycluster%2F",input)

write.table(output.scripts, paste0(path,"/great_jobs.txt"),quote=F,row.names=F,col.names=F)

