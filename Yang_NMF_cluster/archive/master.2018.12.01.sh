mat_path=../../data/snATAC/matrix_Olivier/
out_path=../../analysis/Yang_NMF_method/frontal_cortex/

python snATAC_YLee/snATAC.nmf.lite.py -i $mat_path/tmp.repl1_frontal_cortex_all.sparse.npz \
      -x $mat_path/tmp.repl1_frontal_cortex_all.xgi \
      -y $mat_path/tmp.repl1_frontal_cortex_all.ygi \
      -o $out_path/frontal_cortex \
      -r 10 -n 1 -p 0.05 -c 0 > $out_path/$out_path.frontal_cortex.log


Rscript snATAC_YLee/snATAC.plotH.R -i $out_path/frontal_cortex.H.mx -o $out_path/frontal_cortex.H
Rscript snATAC_YLee/snATAC.plotW.R -i $out_path/frontal_cortex.W.mx -o $out_path/frontal_cortex.W

python snATAC_YLee/snATAC.nmf.stat.py -m $out_path/frontal_cortex.npz \
      -x $out_path/frontal_cortex.xgi \
      -y $out_path/frontal_cortex.ygi \
      --basis $out_path/frontal_cortex.W.mx \
      --coef $out_path/frontal_cortex.H.mx \
      -c 0.2 -o $out_path/frontal_cortex

# calculate cell sparseness and entropy using the statH file
Rscript snATAC_YLee/snATAC.statBox.R -i $out_path/frontal_cortex.statH \
      -o $out_path/frontal_cortex >> $out_path/frontal_cortex.sta.txt

# calculate silhouette and plot tSNE using coefficient matrix H

python snATAC_YLee/snATAC.nmf.plot.py --normH $out_path/frontal_cortex.normH \
      --statH $out_path/frontal_cortex.statH \
      -p 20 -o $out_path/frontal_cortex

# add additional column (age) to the statH file 
Rscript annotate_age_and_cluster_to_cells.r


# make the bam files for each cluster
mkdir $out_path/bams/
python snATAC_YLee/snATAC.nmf.bam.py \
  --bam ../../data/snATAC/bam_bowtie2_Olivier/tmp.repl1_frontal_cortex_all.sorted.bam \
  --statH $out_path/frontal_cortex.statH -o $out_path/bams/frontal_cortex


for i in {1..10}; do 
(
echo $i 
samtools sort $out_path/bams/frontal_cortex.metacell_${i}.bam -o $out_path/bams/frontal_cortex.metacell_${i}.sorted.bam
samtools index $out_path/bams/frontal_cortex.metacell_${i}.sorted.bam 
bamCoverage --bam $out_path/bams/frontal_cortex.metacell_${i}.sorted.bam --outFileFormat bigwig --outFileName $out_path/bams/frontal_cortex.metacell_${i}.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM
) &
done

wait 

cd $out_path/bams && python ~/software/github/seq-min-scripts/make_IGV_session.py mm10 http://renlab.sdsc.edu/yanxiao/mouse_aging/analysis/Yang_NMF_method/frontal_cortex/bams/ test.xml && cd -

## feature counts
awk -v OFS='\t' '{ print $1":"$2"-"$3,$1,$2,$3,"+"}' $mat_path/tmp.repl1_frontal_cortex_all.ygi > $mat_path/tmp.repl1_frontal_cortex_all.saf
featureCounts -a $mat_path/tmp.repl1_frontal_cortex_all.saf -o $out_path/frontal_cortex.counts $out_path/bams/*.sorted.bam -F SAF -T 16

featureCounts -a /mnt/silencer2/home/shz254/annotations/mm10/gencode.vM10.annotation.gene.tss1k.saf -o $out_path/frontal_cortex.prom.counts $out_path/bams/*.sorted.bam -F SAF -T 16


