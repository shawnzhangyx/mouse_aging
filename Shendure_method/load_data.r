setwd("../../analysis/Shendure_cluster/")

library(Matrix)
        library(Matrix)
        library(proxy)
        library(gplots)
        library(Rtsne)
        library(densityClust)
        library(irlba)

tissue="dorsal_hippocampus"
N_READS = 300

my_matrix = readMM(paste0("../../data/snATAC/ct2bin/",tissue,".all.MM.mtx.gz"))

rownames(my_matrix) = read.delim(paste0('../../data/snATAC/ct2bin/',tissue,".all.xgi"), header=FALSE)$V1

tmp = read.delim(paste0('../../data/snATAC/ct2bin/',tissue,".all.ygi"), header=FALSE)
cnames = paste0(tmp$V1,"_",tmp$V2,"_",tmp$V3)
colnames(my_matrix) = cnames
# t transform the matrix
flies_6to8 = t(my_matrix)

#number of fragments per cell. 
num_frag_ncounted = colSums(flies_6to8)
# only keep cells with more than 300 reads. 
flies_6to8 = flies_6to8[,which(num_frag_ncounted > N_READS)]

num_cells_ncounted = rowSums(flies_6to8)
options(repr.plot.width=4, repr.plot.height=4)
hist(log10(num_cells_ncounted),main="No. of Cells Each Site is Observed In",breaks=50)
abline(v=log10(num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[20000]]),lwd=2,col="indianred")


ncounts = flies_6to8[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[20000]]),]
new_counts = colSums(ncounts)
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=0.1)),lwd=2,col="indianred")


# remove the bottom 10% cells. 
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=0.1)]
ncounts = ncounts[rowSums(ncounts) > 0,]

# transform into TF-IDF
nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
idf = as(log(1 + ncol(ncounts) / Matrix::rowSums(ncounts)), "sparseVector")
tf_idf_counts = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs


    #This step can take a minute
    set.seed(0) #For reproducibility
    SVD = irlba(tf_idf_counts, 6, 6)
    sk_diag = matrix(0, nrow=6, ncol=6)
    diag(sk_diag) = SVD$d
    sk_diag[1,1] = 0
    
    LSI_out = t(t(sk_diag %*% t(SVD$v)) %*% t(SVD$u))
    LSI_out = t(scale(t(LSI_out)))
    LSI_out[LSI_out > 1.5] = 1.5
    LSI_out[LSI_out < -1.5] = -1.5


