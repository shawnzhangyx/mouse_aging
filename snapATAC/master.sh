#tissue=DH

## preprocess the bam files into snap PRE format. 
bash -x run_snakemake_snap_pre.sh

## merge each sample in one tissues. 
## perform clustering. 
# manual 
bash -x merge_and_cluster_one_tissue.sh 

## generate the bam files. 
mkdir ../../analysis/snapATAC/$tissue/bam.cluster_age_rep
python split_bam_files.py \
  --tissue $tissue \
  --bam-prefix ../../data/snATAC/bam.filter/  \
  --bam-suffix .filter.bam \
  --statH ../../analysis/snapATAC/$tissue/${tissue}.pool.barcode.meta_info.txt \
  -o ../../analysis/snapATAC/$tissue/bam.cluster_age_rep/${tissue} \
  -p 30

## split bam files
bash -x ./run_snakemake_gen_bam_bigWig.sh

## subsample bam files 
# bash -x subsample_bam_files.sh

## plot cluster specific QC stats.
for tissue in DH HT LM FC BM; do
Rscript plot_celltype_fraction.r ../../analysis/snapATAC/$tissue/${tissue}.pool.barcode.meta_info.txt ../../analysis/snapATAC/$tissue/${tissue}.celltype.fraction.pdf
Rscript plot_cluster_specific_qc.r $tissue
Rscript plot_umap_qc.r $tissue
done

## get the rank for each tissue. 
function get_rank (){
  local tissue=$1
  echo $(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pool.barcode.meta_info.txt|cut -f 9|sort -k1,1nr|head -n1)
  }
echo $(get_rank DH)

## get the cluster specific genes & peaks. 
for tissue in HT DH LM FC BM; do
  rank=$(get_rank $tissue)
  echo $tissue $rank
  Rscript output_cluster_feature.edger.v2.r $tissue $rank
  done
#### special script. Identify the difference between any two clusters. 
# Rscript SPL.compare_cluster_feature.edger.r 

## differential analysis.
#for tissue in DH HT LM FC; do
#  rank=$(get_rank $tissue)
#  echo $tissue $rank
#  Rscript cluster.difftest.r $tissue $rank
#  done
# Differential peak analysis using edgeR. 
Rscript cluster.difftest.edgeR.r DH DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cluster.difftest.edgeR.r FC FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cluster.difftest.edgeR.r HT HT.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cluster.difftest.edgeR.r LM LM.pool.snapATAC.Frag500.TSS7.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
Rscript cluster.difftest.edgeR.r BM BM.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.cluster.RData

## combine all up and down regulated peaks during aging. 
cat ../../analysis/snapATAC/*/age_diff_edgeR.snap/*.up.bed > ../../analysis/snapATAC/all_celltypes/allcelltypes.age_up.peaks.bed
cat ../../analysis/snapATAC/*/age_diff_edgeR.snap/*.down.bed > ../../analysis/snapATAC/all_celltypes/allcelltypes.age_down.peaks.bed




## see if there is any correlation between number of diff cluster and number of cells. 
#Rscript plot_num_diff_peak_vs_num_cells_per_cluster.r 

## find motifs in each set of differential peaks 
for tissue in FC HT LM BM; do 
cd ../../analysis/snapATAC/$tissue/age_diff_edgeR.snap/
mkdir -p motif.homer motif.homer.bg
for file in *{up,down}.bed; do 
  findMotifsGenome.pl $file mm10 motif.homer/${file/.bed/.homer} -nomotif &
done 
wait

## use all peaks as the background. 
for file in *{up,down}.bed; do
findMotifsGenome.pl $file mm10 motif.homer.bg/${file/.bed/.homer} -nomotif -bg ../../../../data/snATAC/peaks/${tissue}_summits.ext1k.bed &
done
wait 
cd -
Rscript combine_de_peak_motif.results.r $tissue
done
## merge the table of differential motif results. 

### GREAT pathway enrichment analysis. 
for tissue in DH FC HT LM BM; do 
echo $tissue 
Rscript great.make_job_list.r $tissue
cd ../../analysis/snapATAC/$tissue/age_diff_edgeR.snap/
~/software/great/greatBatchQuery.py great_jobs.txt
for file in $(ls great_chipseq/*.great.tsv); do 
cut -f 1-22 $file > ${file/.tsv/.red.tsv} 
done
cd -
Rscript great.process_result.r $tissue
done 

########## All celltype processing #######
mkdir ../../analysis/snapATAC/all_celltypes
cd ../../analysis/snapATAC/
bams=*/bam.cluster/*.sorted.bam
featureCounts -a ../../data/snATAC/peaks/all_tissue.merged.peaks.saf -o all_celltypes/all_celltypes.counts $bams -F SAF -T 20 -O
cd -
Rscript all_celltype.clustering.r
