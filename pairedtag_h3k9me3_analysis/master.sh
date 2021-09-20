#merge bams. 
bash merge_bam.sh
bash merge_bam.rep.sh

## prepare bed files for SICER peak calling. 
bamToBed -i bam_split/PT_Ageing_DNA_merge_sorted.10.bam > bed/PT_Ageing_DNA_merge_sorted.10.bed 
# call peak
sh ~/software/SICER/SICER-rb.sh ./bed/ PT_Ageing_DNA_merge_sorted.10.bed ./peaks mm10 1 1000 300 0.8 3000 0.01
sh ~/software/SICER/SICER-rb.sh ./bed/ PT_Ageing_DNA_merge_sorted.10.bed ./peaks mm10 1 5000 300 0.8 10000 0.01

# generate bigWig
bamCoverage -p 30 -e 100 --binSize 25 -b PT_Ageing_DNA_merge_sorted.10.bam -o PT_Ageing_DNA_merge.25.bw --normalizeUsing RPKM


#featureCounts
cd ../../analysis/paired_tag_h3k9me3 
awk -v OFS="\t" '{print "peak"NR,$1,$2,$3,"."}' peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.bed > peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf

featureCounts -a peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf -o SICER-W5000-G10000.read.counts $(ls bam_merge_rep/*.bam) -F SAF -T 16

featureCounts -a peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf -o SICER-W5000-G10000.reps.read.counts $(ls bam_merge/*.bam) -F SAF -T 16

#run differential analysis. 
Rscript process_h3k9me3domain.edger.r

# count RNA reads into h3k9me3 domains and run differential analysis. 
featureCounts -a peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf -o SICER-W5000-G10000.FC.RNA.read.counts $(ls /projects/ps-renlab/lamaral/projects/aging_RNA/FC/data/snRNA_deep_2/split_bams/consistent_clusters/*.bam) -F SAF -T 16
featureCounts -a peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf -o SICER-W5000-G10000.DH.RNA.read.counts $(ls /projects/ps-renlab/lamaral/projects/aging_RNA/DH/data/split_bams/consistent_clusters/*.bam) -F SAF -T 16 
# differential analysis of RNAs in these regions. 
Rscript process_rna_at_h3k9me3domain.edger.r


## test if the down-regulation of H3K9me3 domains are happening in every cells or just in a subset of cells. 

# split the bam into individual cells. 
#for file in bam_merge_rep/L23.03.bam bam_merge_rep/L23.18.bam; do 
ulimit -n 2048
for sample in L23.03.rep1 L23.03.rep2 L23.18.rep1 L23.18.rep2; do 
  echo $sample
  mkdir -p ../../analysis/paired_tag_h3k9me3/bam_split_cell/$sample
  python split_bam_by_cell.py \
    --file ../../analysis/paired_tag_h3k9me3/bam_merge/$sample.bam \
    --outPrefix ../../analysis/paired_tag_h3k9me3/bam_split_cell/$sample/
done

cd ../../analysis/paired_tag_h3k9me3/ 
for sample in L23.03.rep1 L23.03.rep2 L23.18.rep1 L23.18.rep2; do
  files=$(ls bam_split_cell/$sample/*.bam)
#  featureCounts -a peaks/PT_Ageing_DNA_merge_sorted.10-W5000-G10000-E0.01.scoreisland.saf -o bam_split_cell/$sample.read.counts  $files -F SAF -T 8 &
  featureCounts -a h3k9me3_changed.regions.L23.manual.saf -o bam_split_cell/$sample.L23_changed_domain.counts $files -F SAF -T 4 & 
  done

wait


