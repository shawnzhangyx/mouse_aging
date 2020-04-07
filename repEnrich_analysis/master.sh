# split the meta info
tissue=FC
mkdir ../../analysis/repeats_RepEnrich2/metaInfo/ ../../analysis/repeats_RepEnrich2/fastq.cluster
for tissue in DH HT LM; do 
for sample in 03_rep1 03_rep2 10_rep1 10_rep2 18_rep1 18_rep2; do
  grep ${tissue}_${sample} ../../analysis/snapATAC/$tissue/$tissue.pool.barcode.meta_info.txt > ../../analysis/repeats_RepEnrich2/metaInfo/${tissue}_$sample.meta.txt
  done
done

# extract fastq file based on cluster info.
#for sample in 03_rep1 03_rep2 10_rep1 10_rep2 18_rep1 18_rep2; do
#  for rd in R1 R2; do
#  python ~/software/github/single-cell-utils/split_fastq_by_cluster.py \
#  --meta ../../analysis/repeats_RepEnrich2/metaInfo/${tissue}_${sample}.meta.txt \
#  -i ../../data/snATAC/fastq.demultiplex/${tissue}_${sample}/${tissue}_${sample}.${rd}.fastq.gz \
#  --barcode-pos 2 --cluster-pos 9 --fastq-barcode-re "@(.*?):.*" \
#  --outPrefix ../../analysis/repeats_RepEnrich2/fastq.cluster/${tissue}_${sample} \
#  --outAffix ${rd}.fastq.gz &
#  done
#  done


# run snakemake 
snakemake --snakefile Snakemake.RepEnrich2 -j 200 --ri -k \
  --cluster "qsub -l nodes=1:ppn={threads} -N {rule} -q hotel" \
  --jobscript ./jobscript.pbs --jobname "{rulename}.{jobid}.pbs"

