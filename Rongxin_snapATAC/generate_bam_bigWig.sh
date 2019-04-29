tissue=$1
rank=$2

out_path=../../analysis/Rongxin_snapATAC/$tissue/

mkdir -p $out_path/bams $out_path/bigWig.cluster_age_rep $out_path/bigWig.cluster

python split_bam_files.py \
  --bam ../../../data/snATAC/bam_bowtie2_Olivier/filter_bam/${tissue} \
  --statH $out_path/$tissue.pooled.barcode.cluster.stage.rep.txt
  -o $out_path/bams/${tissue}

## sort and generate bigWig for each file.

for i in $(seq 1 $rank); do
    (
    for age in 03 10 18; do
      for rep in rep1 rep2; do
    samtools sort -m 4G $out_path/bams/${tissue}.metacell_${i}.$age.$rep.bam -o $out_path/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam
    samtools index $out_path/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam
    bamCoverage --bam $out_path/bams/${tissue}.metacell_${i}.$age.$rep.sorted.bam --outFileFormat bigwig --outFileName $out_path/bigWig.cluster_age_rep/${tissue}.metacell_${i}.$age.$rep.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM -p 1
      done
    done
    ) &
done

wait

for i in $(seq 1 $rank); do
    (
    echo $i
    samtools merge -f $out_path/bams/${tissue}.metacell_${i}.bam $out_path/bams/${tissue}.metacell_${i}.??.rep?.bam

    samtools sort -m 4G $out_path/bams/${tissue}.metacell_${i}.bam -o $out_path/bams/${tissue}.metacell_${i}.sorted.bam
    samtools index $out_path/bams/${tissue}.metacell_${i}.sorted.bam
    ) &
    done

wait

for i in $(seq 1 $rank); do
    bamCoverage --bam $out_path/bams/${tissue}.metacell_${i}.sorted.bam --outFileFormat bigwig --outFileName $out_path/bigWig.cluster/${tissue}.metacell_${i}.sorted.rpkm.bw --binSize 25 --normalizeUsingRPKM
    done


    cd $out_path/bigWig.cluster && python ~/software/github/seq-min-scripts/make_IGV_session.py mm10 http://renlab.sdsc.edu/yanxiao/mouse_aging/analysis/Rongxin_snapATAC/${tissue}/bigWig.cluster/ test.xml && cd -

