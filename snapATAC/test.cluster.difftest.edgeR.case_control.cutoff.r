#tissue="DH"
tissue=commandArgs(trailing=T)[1] #"DH"
max_cluster=commandArgs(trailing=T)[2] #cl = 1


setwd(paste0("../../analysis/snapATAC/",tissue,"/age_diff_edgeR.case_control"))
cl =1

for (cl in 1:max_cluster) {

a=fread(paste0(cl,".edger.csv"))

print(c(tissue,cl,length(which(a$efdra<0.2)) ))

}



