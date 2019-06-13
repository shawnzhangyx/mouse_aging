#tissue=DH

## preprocess the bam files into snap PRE format. 
bash -x run_snakemake_snap_pre.sh

## merge each sample in one tissues. 
## perform clustering. 
for tissue in HT LM FC; do 
bash -x merge_and_cluster_one_tissue.sh $tissue &
done 

## obtain rank 
rank=$(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pooled.barcode.cluster.stage.rep.txt|cut -f 2|sort -k1,1nr|head -n1)

## generate the bam files. 
mkdir ../../analysis/snapATAC/$tissue/bam.cluster_age_rep
python split_bam_files.py \
  --tissue $tissue \
  --bam-prefix ../../data/snATAC/bam.filter/  \
  --bam-suffix .filter.bam \
  --statH ../../analysis/snapATAC/$tissue/${tissue}.pooled.barcode.cluster.stage.rep.txt \
  -o ../../analysis/snapATAC/$tissue/bam.cluster_age_rep/${tissue} \
  -p 40

## split bam files
#bash -x generate_bam_bigWig_rankR.sh $tissue $rank
bash -x ./run_snakemake_gen_bam_bigWig.sh

function get_rank (){
  local tissue=$1
  echo $(tail -n+2 ../../analysis/snapATAC/${tissue}/${tissue}.pooled.barcode.cluster.stage.rep.txt|cut -f 2|sort -k1,1nr|head -n1)
  }
echo $(get_rank DH)


## differential analysis.
for tissue in DH HT LM FC; do
  rank=$(get_rank $tissue)
  echo $tissue $rank
  Rscript cluster.difftest.r $tissue $rank
  done
