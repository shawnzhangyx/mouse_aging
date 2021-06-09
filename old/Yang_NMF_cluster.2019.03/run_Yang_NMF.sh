tissue=heart
R=10
mat_path=../../data/snATAC/ct2peaks/
out_path=../../analysis/Yang_NMF_method/${tissue}/
mkdir $out_path

python snATAC_YLee/snATAC.nmf.lite.py -i $mat_path/${tissue}.all.npz \
      -x $mat_path/${tissue}.all.xgi \
      -y $mat_path/${tissue}.all.ygi \
      -o $out_path/${tissue} \
      -r 5 -n 1 -p 0.05 -c 500 > $out_path/${tissue}.log

# plot H & W. 
#Rscript snATAC_YLee/snATAC.plotH.R -i $out_path/${tissue}.H.mx -o $out_path/${tissue}.H
#Rscript snATAC_YLee/snATAC.plotW.R -i $out_path/${tissue}.W.mx -o $out_path/${tissue}.W

python snATAC_YLee/snATAC.nmf.stat.py -m $out_path/${tissue}.npz \
      -x $out_path/${tissue}.xgi \
      -y $out_path/${tissue}.ygi \
      --basis $out_path/${tissue}.W.mx \
      --coef $out_path/${tissue}.H.mx \
      -c 0.2 -o $out_path/${tissue}

# calculate cell sparseness and entropy using the statH file
#Rscript snATAC_YLee/snATAC.statBox.R -i $out_path/${tissue}.statH \
#      -o $out_path/${tissue} >> $out_path/${tissue}.sta.txt

# calculate silhouette and plot tSNE using coefficient matrix H
python snATAC_YLee/snATAC.nmf.plot.py --normH $out_path/${tissue}.normH \
      --statH $out_path/${tissue}.statH \
      -p 200 -o $out_path/${tissue}

# add additional column (age) to the statH file 
Rscript annotate_age_and_cluster_to_cells.r $tissue
# plot age replicate 
Rscript plot_celltype_fraction.r $tissue

## split cells by age and cluster. 
mkdir $out_path/bams
python split_bam_files.py \
  --bam ../../data/snATAC/bam_bowtie2_Olivier/filter_bam/${tissue} \
  --statH $out_path/${tissue}.barcode.cluster.stage.rep.txt -o $out_path/bams/${tissue}


for i in {1..5}; do 
(
echo $i 
samtools merge -f $out_path/bams/${tissue}.metacell_${i}.bam $out_path/bams/${tissue}.metacell_${i}.??.rep?.bam 

samtools sort -m 4G $out_path/bams/${tissue}.metacell_${i}.bam -o $out_path/bams/${tissue}.metacell_${i}.sorted.bam
samtools index $out_path/bams/${tissue}.metacell_${i}.sorted.bam 
bamCoverage --bam $out_path/bams/${tissue}.metacell_${i}.sorted.bam --outFileFormat bigwig --outFileName $out_path/bams/${tissue}.metacell_${i}.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM
) &
done

wait 

cd $out_path/bams && python ~/software/github/seq-min-scripts/make_IGV_session.py mm10 http://renlab.sdsc.edu/yanxiao/mouse_aging/analysis/Yang_NMF_method/${tissue}/bams/ test.xml && cd -

## feature counts
echo "counting the peak regions" 
#awk -v OFS='\t' '{ print $1":"$2"-"$3,$1,$2,$3,"+"}' $mat_path/tmp.repl1_${tissue}_all.ygi > $mat_path/tmp.repl1_${tissue}_all.saf
featureCounts -a ../../data/snATAC/peaks/${tissue}.all.pooled_summits.ext1k.filter.saf -o $out_path/${tissue}.counts $out_path/bams/*.rep?.bam -F SAF -T 16
echo "Counting the promoter regions"
featureCounts -a /mnt/silencer2/home/shz254/annotations/mm10/gencode.vM10.annotation.gene.tss1k.saf -o $out_path/${tissue}.prom.counts $out_path/bams/*.rep?.bam -F SAF -T 16


