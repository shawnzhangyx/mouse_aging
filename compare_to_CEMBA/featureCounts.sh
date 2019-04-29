bams=/projects/ps-renlab/yangli/projects/CEMBA/01.nmf/CEMBA1*[^k]/*r15*/*.bam
out_path=../../analysis/compare_to_CEMBA
mkdir $out_path

featureCounts -a /mnt/silencer2/home/shz254/annotations/mm10/gencode.vM10.annotation.gene.tss1k.saf -o $out_path/CEMBA.prom.counts $bams -F SAF -T 16

