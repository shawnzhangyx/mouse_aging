tissue=heart
tissue=$1
path=../../data/snATAC/bam_bowtie2_Olivier/
mkdir $path/filter_bam/

for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
  samtools view -h $path/dedup_bam/$tissue.$sample.bc.dedup.bam |\
  awk -v OFS="\t" ' function abs(v) {return v < 0 ? -v : v} 
  { if ( /^@/) { print $0} else if ( $7=="=" && abs($9)< 2000) {print $0}}'|\
  samtools view -bS > $path/filter_bam/$tissue.$sample.dedup.filter.bam &
  done

wait

for sample in 03.rep1 03.rep2 10.rep1 10.rep2 18.rep1 18.rep2; do
  samtools view $path/filter_bam/$tissue.$sample.dedup.filter.bam |\
    awk -v OFS="\t" '{split($1,a,":"); print a[1],$9}' |\
    gzip >  $path/insert_size/$tissue.$sample.dedup.filter.insert_size.txt.gz &
done
