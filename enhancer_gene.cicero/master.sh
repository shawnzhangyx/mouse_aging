mkdir ../../analysis/cicero_results/
# running cicero to detect links between peaks. 
Rscript cicero.r DH ../../snapATAC/DH/snapFiles/DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cicero.r FC ../../snapATAC/FC/snapFiles/FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cicero.r HT ../../snapATAC/HT/snapFiles/HT.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cicero.r LM ../../snapATAC/LM/snapFiles/LM.pool.snapATAC.Frag500.TSS7.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
# apply a cutoff on the cicero connections. 
cd ../../analysis/cicero_results/
for tissue in DH FC HT LM; do
grep -v NA ${tissue}/${tissue}.cicero_conns.csv |awk -v FS=',' '{ if (NR==1){print $0} else { if ($4>0.25){print $0}}}'  > ${tissue}/${tissue}.cicero_conns.25.csv
Rscript ../../scripts/enhancer_gene.cicero/convert.cicero_csv_to_bedped.r ${tissue}/${tissue}.cicero_conns.25.csv ${tissue}/${tissue}.cicero_conns.25.bedpe
done
cd -

# overlap peaks that anchored on promoter of longevity genes. 
cd ../../analysis/cicero_results/
for tissue in DH FC HT LM; do 
intersectBed -a ../../data/snATAC/peaks/${tissue}_summits.ext1k.bed -b ../../annotations/gencode.vM10.annotation.gene.tss1k.bed -wo > ${tissue}/${tissue}_peak_overlap_promoter.txt
done
cd -

# find the co-accessible peaks. 
for tissue in DH FC HT LM; do
Rscript find_co_access_peaks.r $tissue
done

# overlap the co-accessible peaks with GWAS SNPs. 
cd ../../analysis/cicero_results/
for tissue in DH FC HT LM; do
intersectBed -a ${tissue}/${tissue}_peak_coaccess_toGene.25.bed -b /projects/ps-renlab/lamaral/projects/Aging/GWAS/bed_LD_snp_files_wtag/hglft_genome_wtag.bed -wo > ${tissue}/${tissue}_coA_peak.overlapGWAS.25.txt 
done

 > all_tissue.coA_peak.overlapGWAS.25.txt
for tissue in DH FC HT LM; do
awk -v OFS="\t" -v tissue=$tissue '{print $0,tissue}'  ${tissue}/${tissue}_coA_peak.overlapGWAS.25.txt >> all_tissue.coA_peak.overlapGWAS.25.txt 
done

