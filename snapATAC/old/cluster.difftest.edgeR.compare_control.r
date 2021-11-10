tissue=commandArgs(trailing=T)[1] #"DH"
max_cluster=commandArgs(trailing=T)[2] #cl = 1

setwd(paste0("../../analysis/snapATAC/",tissue))

library(doParallel)
registerDoParallel(cores=max_cluster)

foreach(cl=1:max_cluster) %dopar% {

test = read.csv(paste0("age_diff_edgeR.snap/",cl,".edger.txt"))
test$F = 0
control = read.csv(paste0("age_diff_edgeR.control/",cl,".edger.txt"))
control$F = 1
com = rbind(test,control)

com = com[order(com$PValue),]

com$rank = 1:nrow(com)
com$F_cum = com$F
for (i in 2:nrow(com)) {
  com$F_cum[i] = com$F_cum[i-1] + com$F[i]
  }


com$efdr = com$F_cum/(com$rank-com$F_cum)
com$efdra = com$efdr
## 
for (i in (nrow(com)-1):1) {
  if ( com$efdra[i] > com$efdra[i+1]) {
  com$efdra[i] = com$efdra[i+1] 
    }
}
write.csv(com,paste0("age_diff_edgeR.case_control/",cl,".edger.csv"))
}

