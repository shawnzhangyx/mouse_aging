link=data.frame(fread("/mnt/tscc/lamaral/projects/aging_hippocampus_RNA/analysis/weighted_regression_500kb+500kb-.txt"))
colnames(link) = c("Geneid","Chr","Start","End","Peakid","ATAC_Chr","ATAC_start","ATAC_end","slope","adj.r.squared","Pvalue")

link$fdr = p.adjust(link$Pvalue,method="BH")

sub = subset(link, Geneid %in% c("Robo1","Apoe","Nrg1") & fdr<0.05)

beds = sub[,c(6,7,8,5,12)]

write.table(beds,"../../analysis/rna_atac_integration/peaks_linked_toGOI.bed",row.names=F,col.names=F,quote=F,sep="\t")

