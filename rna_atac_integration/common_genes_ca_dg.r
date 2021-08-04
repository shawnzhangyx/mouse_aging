setwd("../../analysis/rna_atac_integration/RNA_diff")

ca1 = data.frame(fread("CA1_03vs18_nocutoff.txt"))
ca1 = ca1[which(ca1$p_val_adj < 0.05),]
ca2 = data.frame(fread("CA23_03vs18_nocutoff.txt"))
ca2 = ca2[which(ca2$p_val_adj < 0.05),]

dg = data.frame(fread("DG_03vs18_nocutoff.txt"))
dg = dg[which(dg$p_val_adj < 0.05),]

se = data.frame(fread("Sub_Ent_03vs18_nocutoff.txt"))
se = se[which(se$p_val_adj < 0.05),]

> table(ca1$V1 %in% ca2$V1)

FALSE  TRUE
  233    21

> table(ca1$V1 %in% dg$V1)
FALSE  TRUE
  171    83


com = ca1$V1[which(ca1$V1 %in% dg$V1)]

ca2 = data.frame(fread("CA23_03vs18_nocutoff.txt"))

out = data.frame(com,ca1[match(com,ca1$V1),c(3,6)],
    ca2[match(com,ca2$V1),c(3,6)],
    dg[match(com,dg$V1),c(3,6)],
    se[match(com,se$V1),c(3,6)]
    )

out[which(out$p_val_adj.1<0.05),]





