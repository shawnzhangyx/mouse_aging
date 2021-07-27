a=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/final_pairs_table_concat.txt",skip=1))



setwd("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/snapATAC/DH/age_diff_edgeR.snap/")

files = list.files(pattern="*(up|down).bed")

dat.list = list()

for (file in files) { 
  dat.list[[file]] = data.frame(fread(file))
  dat.list[[file]]$file = file
  }

dat = do.call(rbind,dat.list)

dat$peak = paste(dat$V1,dat$V2,dat$V3)

a$peak = paste(a$V1,a$V2,a$V3)



## overlap with diff peaks. 
b = a[which(a$peak %in% dat$peak),]
peaks = dat[which(dat$peak %in% a$peak),]
peaks$gene = b$V8[match(peaks$peak,b$peak)]
#                                  peak gene
#7.up.bed.1294  chr7 19823890 19824891 Apoe
#8.down.bed.212 chr7 19875634 19876635 Apoe
#8.down.bed.213 chr7 19876025 19877026 Apoe
#9.down.bed.48  chr7 19600165 19601166 Apoe

setwd("/projects/ps-renlab/lamaral/projects/aging_hippocampus_RNA/analysis/diff_exp_after_doub_rem/")


fs2 = list.files(pattern="3vs18.txt")

dat.list = list()

for (file in fs2) {
  dat.list[[file]] = data.frame(fread(file))
  dat.list[[file]]$file = file
    }
dat = do.call(rbind,dat.list)

dat2 = dat[which(dat$p_val_adj<0.05),]

dat2$p_val = -log10(dat2$p_val)
dat2$p_val_adj = -log10(dat2$p_val_adj)

genes = unique(dat2$V1)

dat2 = dat2[which(dat2$V1 %in% b$V8),]

# Promoter and gene both change. 
dat2[which(dat2$V1=="Lin28b"),]
peaks[which(peaks$gene=="Lin28b"),]


# potential enhancer at Robo1 in DG. 
dat2[which(dat2$V1=="Robo1"),]
peaks[which(peaks$gene=="Robo1"),]


Pcdh7
dat2[which(dat2$V1=="Nrxn3"),]
peaks[which(peaks$gene=="Nrxn3"),]

