#tissue=heart
#rank=10
tissue=$1
rank=$2

out_path=../../analysis/Yang_NMF_method/${tissue}/

mkdir $out_path/R$rank/bams $out_path/R$rank/bigWig.cluster $out_path/R$rank/bigWig.cluster_age bigWig.cluster_age_rep

python split_bam_files.py \
  --bam ../../data/snATAC/bam_bowtie2_Olivier/filter_bam/${tissue} \
    --statH $out_path/R$rank/${tissue}.R$rank.barcode.cluster.stage.rep.txt -o $out_path/R$rank/bams/${tissue}

## sort and generate bigWig for each file. 

for i in $(seq 1 $rank); do
    (
    for age in 03 10 18; do
      for rep in rep1 rep2; do 
    samtools sort -m 4G $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.$rep.bam -o $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam
    samtools index $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam
    bamCoverage --bam $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam --outFileFormat bigwig --outFileName $out_path/R$rank/bigWig.cluster_age_rep/${tissue}.metacell_${i}.$age.$rep.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM -p 4 
      done
    done
    ) & 
done


# merge the files in the same cluster
for i in $(seq 1 $rank); do
    (
    echo $i
    samtools merge -f $out_path/R$rank/bams/${tissue}.metacell_${i}.bam $out_path/R$rank/bams/${tissue}.metacell_${i}.??.rep?.bam

    samtools sort -m 4G $out_path/R$rank/bams/${tissue}.metacell_${i}.bam -o $out_path/R$rank/bams/${tissue}.metacell_${i}.sorted.bam
    samtools index $out_path/R$rank/bams/${tissue}.metacell_${i}.sorted.bam
    ) &
    done

wait

for i in $(seq 1 $rank); do
    bamCoverage --bam $out_path/R$rank/bams/${tissue}.metacell_${i}.sorted.bam --outFileFormat bigwig --outFileName $out_path/R$rank/bigWig.cluster/${tissue}.metacell_${i}.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM
    done



# merge the files in the same cluster and age

for i in $(seq 1 $rank); do
    (
  for age in 03 10 18; do 
    echo $i
    samtools merge -f $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.bam $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.rep?.bam

    samtools sort -m 4G $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.bam -o $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.sorted.bam
    samtools index $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.sorted.bam
  done
    ) &
    done

wait

for i in $(seq 1 $rank); do
    for age in 03 10 18; do
    bamCoverage --bam $out_path/R$rank/bams/${tissue}.metacell_${i}.$age.sorted.bam --outFileFormat bigwig --outFileName $out_path/R$rank/bigWig.cluster_age/${tissue}.metacell_${i}.$age.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM
    done
done



    cd $out_path/R$rank/bigWig.cluster && python ~/software/github/seq-min-scripts/make_IGV_session.py mm10 http://renlab.sdsc.edu/yanxiao/mouse_aging/analysis/Yang_NMF_method/${tissue}/R$rank/bigWig.cluster/ test.xml && cd -
    cd $out_path/R$rank/bigWig.cluster_age && python ~/software/github/seq-min-scripts/make_IGV_session.py mm10 http://renlab.sdsc.edu/yanxiao/mouse_aging/analysis/Yang_NMF_method/${tissue}/R$rank/bigWig.cluster_age/ test.xml && cd -


    ## feature counts
    echo "counting the peak regions"
    #awk -v OFS='\t' '{ print $1":"$2"-"$3,$1,$2,$3,"+"}' $mat_path/tmp.repl1_${tissue}_all.ygi > $mat_path/tmp.repl1_${tissue}_all.saf
    featureCounts -a ../../data/snATAC/peaks/${tissue}.all.pooled_summits.ext1k.filter.saf -o $out_path/R$rank/${tissue}.counts $out_path/R$rank/bams/*.rep?.bam -F SAF -T 16
    echo "Counting the promoter regions"
    featureCounts -a /mnt/silencer2/home/shz254/annotations/mm10/gencode.vM10.annotation.gene.tss1k.saf -o $out_path/R$rank/${tissue}.prom.counts $out_path/R$rank/bams/*.rep?.bam -F SAF -T 16



