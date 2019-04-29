#!/bin/sh                                                                                                                                                                                                

dropseq_dir=/mnt/silencer1/home/amraman/Sebastian/Dropseq_intron
picard_dir=/mnt/silencer1/home/amraman/Sebastian/Dropseq_intron/jar
star_dir=/mnt/silencer1/home/amraman/STAR-2.5.2b/bin/Linux_x86_64/
ref_dir=/mnt/silencer1/home/amraman/Sebastian/mouse_gencode
tmp_path=`pwd`
#############################################
sampleName=JB_55
species='m'
exonOnly='F'
cells=5000 
############################################

while getopts ":n:s:e:c:" options; do
    case $options in
        n ) sampleName=$OPTARG;;
        s ) species=$OPTARG;;
        e ) exonOnly=$OPTARG;;
        c ) cells=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

# combine lane 1 and lane 2 fastq files together                                                                                                                                                                                                                                                                                             
echo "Start to cat files ..."                                                                                                                                                                          
##########################################
####### if files are in pwd ##############
cat *R1* > $sampleName.r1.fastq.gz
cat *R2* > $sampleName.r2.fastq.gz
#############################################

read1=$sampleName.r1.fastq.gz 
read2=$sampleName.r2.fastq.gz 
mkdir -p Reports
mkdir -p Tmp
############################################

#extract cell barcodes and UMI from read1, convert pair-end fastq files to bam file, and generate fastq file for alignment                                                                               
java -Xmx4g -jar $picard_dir/picard.jar FastqToSam FASTQ=$read1 FASTQ2=$read2 SAMPLE_NAME=$sampleName OUTPUT=/dev/stdout TMP_DIR='pqd'/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=$sampleName'.unaligned.sorted.bam' SORT_ORDER=queryname TMP_DIR=`pwd`/Tmp 

######cell barcode (bases1-16); umi barcode (17-26) 
$dropseq_dir/TagBamWithReadSequenceExtended I=$sampleName'.unaligned.sorted.bam' O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.cell_tag_report.txt' BASE_RANGE=1-16 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=`pwd`/Reports/$sampleName'.molecule_tag_report.txt' BASE_RANGE=17-26 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/FilterBAM TAG_REJECT=XQ I=/dev/stdin O=/dev/stdout | \
$dropseq_dir/TrimStartingSequence INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=`pwd`/Reports/$sampleName'.adapter_trimming_report.txt' SEQUENCE='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' MISMATCHES=0 NUM_BASES=5 | \
$dropseq_dir/PolyATrimmer INPUT=/dev/stdin OUTPUT=/dev/stdout MISMATCHES=0 NUM_BASES=6 | \
tee $sampleName'.unaligned.tagged.bam' | \
java -Xmx8g -jar $picard_dir/picard.jar SamToFastq INPUT=/dev/stdin FASTQ=`pwd`/$sampleName'.unaligned.tagged.fastq'

###
#select ref index genome generate by STAR, and perform STAR alignment                                                                                                                                     
if [ $species = 'h' ]; then
        refSTAR=humGenomeIndex_99
        refFasta=GRCh38.primary_assembly.genome.fa 
                refGTF=gencode.v24.primary_assembly.annotation.gtf
else
        refSTAR=mouse_gencode
        refFasta=GRCm38.primary_assembly.genome.fa
                refGTF=gencode.vM15.annotation.gtf
fi




$star_dir/STAR --runThreadN 16 --genomeDir $ref_dir --readFilesIn `pwd`/$sampleName'.unaligned.tagged.fastq' --outFileNamePrefix $sampleName --outReadsUnmapped Fastx
#merge aligned sam file with cell barcode/UMI tagged bam file, correct barcode synthesis error, and generate digital expression matrix                                                                   
mkdir -p DGE
#############
## take the mapped reads and map them to exon                                                                                                                                                 
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=$sampleName'.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=$tmp_path/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout PAIRED_RUN=false |\
$dropseq_dir/TagReadWithGeneExon I=/dev/stdin O=$sampleName'.aligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE

## take the unmapped reads and map them to gene                                                                                                                                               
$star_dir/STAR --runThreadN 16 --genomeDir $ref_dir --readFilesIn `pwd`/$sampleName'.Unmapped.out.mate1' --outFileNamePrefix $sampleName'.map2' --outSAMunmapped Within                                  
java -Xmx4g -jar $picard_dir/picard.jar SortSam I=$sampleName'.map2.Aligned.out.sam' O=/dev/stdout SO=queryname TMP_DIR=$tmp_path/Tmp |\
java -Xmx4g -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout PAIRED_RUN=false |\
$dropseq_dir/TagReadWithGene I=/dev/stdin O=$sampleName'.intronaligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE

java -Xmx4g -jar $picard_dir/picard.jar MergeSamFiles I=$sampleName'.intronaligned.gene.bam' I=$sampleName'.aligned.gene.bam' O=$sampleName'.aligned_both.gene.bam'
$dropseq_dir/BAMTagHistogram I=$sampleName'.aligned_both.gene.bam' O=`pwd`/DGE/$sampleName'.intron_fixreadsByBarcode.txt.gz' TAG=XC
$dropseq_dir/DigitalExpression I=$sampleName'.aligned_both.gene.bam' O=`pwd`/DGE/$sampleName'.intron_fix_counts.tsv' SUMMARY=`pwd`/Reports/$sampleName'.intron_fix.count_summary.txt' NUM_CORE_BARCODES=$cells EDIT_DISTANCE=1 

##########################################
mv $sampleName'.Log.out' Reports/
mv $sampleName'.Log.progress.out' Reports/
mv $sampleName'.SJ.out.tab' Reports/
mv $sampleName'.Log.final.out' Reports/

mv $sampleName'.unaligned.sorted.bam' Tmp/
mv $sampleName'.aligned.gene.bam' Tmp/
mv $sampleName'.aligned.clean.bam' Tmp/
mv $sampleName'.unaligned.tagged.fastq' Tmp/
mv $sampleName'.unaligned.tagged.bam' Tmp/
mv $sampleName'.Aligned.out.sam' Tmp/



