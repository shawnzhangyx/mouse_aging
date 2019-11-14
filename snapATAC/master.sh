#tissue=DH

## preprocess the bam files into snap PRE format. 
bash -x run_snakemake_snap_pre.sh

## merge each sample in one tissues. 
## perform clustering. 
for tissue in HT LM FC; do 
bash -x merge_and_cluster_one_tissue.sh $tissue &
done 

## obtain rank 
rank=$(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pool.barcode.meta_info.txt|cut -f 11|sort -k1,1nr|head -n1)

## generate the bam files. 
mkdir ../../analysis/snapATAC/$tissue/bam.cluster_age_rep
python split_bam_files.py \
  --tissue $tissue \
  --bam-prefix ../../data/snATAC/bam.filter/  \
  --bam-suffix .filter.bam \
  --statH ../../analysis/snapATAC/$tissue/${tissue}.pool.barcode.meta_info.txt \
  -o ../../analysis/snapATAC/$tissue/bam.cluster_age_rep/${tissue} \
  -p 40

## split bam files
#bash -x generate_bam_bigWig_rankR.sh $tissue $rank
bash -x ./run_snakemake_gen_bam_bigWig.sh

## subsample bam files 
bash -x subsample_bam_files.sh


## plot cluster specific QC stats.
for tissue in DH HT LM FC; do
Rscript plot_cluster_specific_qc.r $tissue
done

## get the rank for each tissue. 
function get_rank (){
  local tissue=$1
  echo $(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pooled.barcode.cluster.stage.rep.txt|cut -f 2|sort -k1,1nr|head -n1)
  }
echo $(get_rank DH)

## get the cluster specific genes & peaks. 
for tissue in HT DH LM FC; do
  rank=$(get_rank $tissue)
  echo $tissue $rank
  Rscript output_cluster_feature.edger.r $tissue $rank
  done
#### special script. Identify the difference between any two clusters. 
# Rscript SPL.compare_cluster_feature.edger.r 

## differential analysis.
for tissue in DH HT LM FC; do
  rank=$(get_rank $tissue)
  echo $tissue $rank
  Rscript cluster.difftest.r $tissue $rank
  done

## see if there is any correlation between number of diff cluster and number of cells. 
Rscript plot_num_diff_peak_vs_num_cells_per_cluster.r 






