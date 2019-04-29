cd ../../data/snATAC
bowtie2  --very-sensitive  -x /projects/ps-renlab/share/bowtie2_indexes/mm10 \
  -1 demultiplexed_Olivier/frontal_cortex_all.demultiplexed.R1.fastq.gz \
  -2 demultiplexed_Olivier/frontal_cortex_all.demultiplexed.R2.fastq.gz \
  -p 16 -X 2000 \
  |  samtools view -u -  \
  |  samtools sort -  > bowtie2_bam/frontal_cortex_all.sorted.bam

cd ../../data/snATAC
bowtie2  --very-sensitive  -x /projects/ps-renlab/share/bowtie2_indexes/mm10 \
  -1 demultiplexed_Olivier/frontal_cortex_3_months_v2.demultiplexed.R1.fastq.gz \
  -2 demultiplexed_Olivier/frontal_cortex_3_months_v2.demultiplexed.R2.fastq.gz \
  -p 16 -X 2000 \
  |  samtools view -u -  \
  |  samtools sort -  > bowtie2_bam/frontal_cortex_3_months.sorted.bam

cd ../../data/snATAC
bowtie2  --very-sensitive  -x /projects/ps-renlab/share/bowtie2_indexes/mm10 \
  -1 demultiplexed_Olivier/frontal_cortex_10_months.demultiplexed.R1.fastq.gz \
  -2 demultiplexed_Olivier/frontal_cortex_10_months.demultiplexed.R2.fastq.gz \
  -p 16 -X 2000 \
  |  samtools view -u -  \
  |  samtools sort -  > bowtie2_bam/frontal_cortex_10_months.sorted.bam

cd ../../data/snATAC
bowtie2  --very-sensitive  -x /projects/ps-renlab/share/bowtie2_indexes/mm10 \
  -1 demultiplexed_Olivier/frontal_cortex_18_months.demultiplexed.R1.fastq.gz \
  -2 demultiplexed_Olivier/frontal_cortex_18_months.demultiplexed.R2.fastq.gz \
  -p 16 -X 2000 \
  |  samtools view -u -  \
  |  samtools sort -  > bowtie2_bam/frontal_cortex_18_months.sorted.bam

