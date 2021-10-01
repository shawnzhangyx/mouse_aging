#setwd("../../analysis/rna_atac_integration/RNA_diff")
setwd("/projects/ps-renlab/lamaral/projects/aging_RNA/FC/analysis/seurat_DE_results/")
ct1 = c("Asc","Claustrum","Endo","Inh.Lrrc38.Plcxd3","Inh.Meis2","Inh.Sst.Pvalb","Inh.Vip.Npy","L2-3","L4","L5.Deptor","L5.Fezf2","L5.Parm1","L6","Mgc","Ogc","Opc","Peri")
ct2 = c("Glia","ExN","Glia","InN","InN","InN","InN","ExN","ExN","ExN","ExN","ExN","ExN","Glia","Glia","Glia","Glia")

up_list = list()
down_list = list()

for (ct in ct1) {
  tmp = data.frame(fread(paste0(ct,"_03vs18_nocutoff.txt")))
  sig = tmp[which(tmp$p_val_adj <0.05),]
  if (length(which(sig$avg_logFC<0))>0) { 
  up = cbind(sig[which(sig$avg_logFC<0),c(1,3,6)],ct)
  up_list[[ct]] = up 
  }
  if (length(which(sig$avg_logFC>0))>0) {
  down = cbind(sig[which(sig$avg_logFC>0),c(1,3,6)],ct)
  down_list[[ct]] = down }

}

up1 = do.call(rbind,up_list)
down1 = do.call(rbind,down_list)

up1$ct2 = ct2[match(up1$ct,ct1)]
down1$ct2 = ct2[match(down1$ct,ct1)]

## 
setwd("/projects/ps-renlab/yanxiao/projects/mouse_aging/analysis/rna_atac_integration/RNA_diff")

ct1 = c("Asc","CA1","CA23","DG","Endo","Inh","Mgc","Ogc","Opc","Peri","SMC","Sub_Ent")
ct2 = c("Glia","ExN","ExN","ExN","Glia","InN","Glia","Glia","Glia","Glia","Glia","ExN")

up_list = list()
down_list = list()

for (ct in ct1) {
  tmp = data.frame(fread(paste0(ct,"_03vs18_nocutoff.txt")))
  sig = tmp[which(tmp$p_val_adj <0.05),]
  if (length(which(sig$avg_logFC<0))>0) { 
  up = cbind(sig[which(sig$avg_logFC<0),c(1,3,6)],ct)
  up_list[[ct]] = up
  }
  if (length(which(sig$avg_logFC>0))>0) {
  down = cbind(sig[which(sig$avg_logFC>0),c(1,3,6)],ct)
  down_list[[ct]] = down }
  
}

up2 = do.call(rbind,up_list)
down2 = do.call(rbind,down_list)
up2$ct2 = ct2[match(up2$ct,ct1)]
down2$ct2 = ct2[match(down2$ct,ct2)]


## total cell types: 
# ExN Glia  InN 
#  11   13    5 


up = rbind(up1,up2)
down = rbind(down1,down2)
## test the table. 
table(up$ct,up$ct2)
table(down$ct,down$ct2)

up$count = 1
up_count = aggregate(count~V1+ct2,up,sum)
up_count.t = reshape(up_count,direction="wide",idvar="V1",timevar="ct2")
up_count.t[is.na(up_count.t)] = 0
up_count.t = up_count.t[order(-up_count.t$count.ExN),]
up_count.t$ratio.ExN = up_count.t$count.ExN/11
up_count.t$ratio.Glia = up_count.t$count.Glia/13
up_count.t$ratio.InN= up_count.t$count.InN/5


down$count = 1
down_count = aggregate(count~V1+ct2,down,sum)
down_count.t = reshape(down_count,direction="wide",idvar="V1",timevar="ct2")
down_count.t[is.na(down_count.t)] = 0
down_count.t = down_count.t[order(-down_count.t$count.ExN),]
down_count.t$ratio.ExN = down_count.t$count.ExN/11
down_count.t$ratio.Glia = down_count.t$count.Glia/13
down_count.t$ratio.InN= down_count.t$count.InN/5

setwd("../")
write.csv(up_count.t,"genes_up_in_18m.all_celltypes.csv")
write.csv(down_count.t,"genes_down_in_18m.all_celltypes.csv")


