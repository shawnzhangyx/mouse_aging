# Snakefile.pre_process: align fastq to get bam and filter bam files.
snakemake -j30 --snakefile Snakemake.pre_process
#Fix bam mate (for calculating TSS enrichment): 
for sample in FC_30y_rep1 FC_30y_rep2 FC_90y_rep1; do 
  samtools fixmate -m ../../data/human_cortex/bam.filter.nsort/$sample/${sample}.filter.nsort.bam ../../data/human_cortex/TSSvsFrag.out/${sample}.fixmate.bam -@ 8 &
  done

# generate TSS enrichment
for sample in FC_30y_rep1 FC_30y_rep2 FC_90y_rep1; do 
/projects/ps-renlab/yangli/scripts/snATACutils/snapATAC.qc.TSSvsFrag \
  ../../annotations/hg38.gencode.v26.annotation.gtf \
  ../../data/human_cortex/TSSvsFrag.out/${sample}.fixmate.bam \
  ../../data/human_cortex/TSSvsFrag.out/${sample}/ &
  done
# summarize all information for samples. 
for sample in FC_30y_rep1 FC_30y_rep2 FC_90y_rep1; do
  Rscript summalize_bc_info.per_sample.r $sample; 
done
# Snakefile.snap_pre. Run snaptools to get snap file. 
snakemake -j30 --snakefile Snakemake.snap_pre
# merge and cluster (snapATAc)  on each tissue.
bash merge_and_cluster.sh

## generate the split bam files.
mkdir ../../analysis/human_cortex/$tissue/bam.cluster_age_rep
python split_bam_files.py \
  --tissue $tissue \
  --bam-prefix ../../data/human_cortex/bam.filter/  \
  --bam-suffix .filter.bam \
  --statH ../../analysis/human_cortex/$tissue/${tissue}.pool.barcode.meta_info.txt \
  -o ../../analysis/human_cortex/$tissue/bam.cluster_age_rep/${tissue} \
  -p 30


Rscript plot_celltype_fraction.r ../../analysis/human_cortex/$tissue/${tissue}.pool.barcode.meta_info.txt ../../analysis/human_cortex/$tissue/${tissue}.celltype.fraction.pdf
Rscript plot_cluster_specific_qc.r $tissue
Rscript plot_umap_qc.r $tissue


## differential test.
Rscript cluster.difftest.edgeR.r FC FC.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.harmony.cluster.RData
