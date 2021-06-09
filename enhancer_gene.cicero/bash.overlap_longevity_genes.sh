# Curate symbol for longevity genes.
Rscript curate_longevity_genes.r
# find promoter of these genes.
Rscript find_promoter_of_longevity_genes.r
# overlap peaks that anchored on promoter of longevity genes.
cd ../../analysis/cicero_results/
for tissue in DH FC HT LM; do
intersectBed -a ../../data/snATAC/peaks/${tissue}_summits.ext1k.bed -b genage.tss1k.bed -wo > ${tissue}/${tissue}_peak_overlap_genage_promoter.txt
done
cd -

# find the co-accessible peaks.
for tissue in DH FC HT LM; do
Rscript find_co_access_genage_peaks.r $tissue
done

# overlap the co-accessible peaks with GWAS SNPs.
cd ../../analysis/cicero_results/
for tissue in DH FC HT LM; do
intersectBed -a ${tissue}/${tissue}_peak_coaccess_toGene.25.bed -b /projects/ps-renlab/lamaral/projects/Aging/GWAS/hglft_genome_all_LD_snps.bed -wo > ${tissue}/${tissue}_coA_peak.overlapGWAS.25.txt
done

 > all_tissue.coA_peak.overlapGWAS.25.txt
for tissue in DH FC HT LM; do
awk -v OFS="\t" -v tissue=$tissue '{print $0,tissue}'  ${tissue}/${tissue}_coA_peak.overlapGWAS.25.txt >> all_tissue.coA_peak.overlapGWAS.25.txt
done

