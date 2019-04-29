tissue=dorsal_hippocampus

# remove duplicates.
bash -x ./rmdup.sh $tissue

# keep only reads that map within 2000 bp of the mate pair. 
bash -x filter_bam_by_insert_size.sh $tissue

# perform qc. 
bash -x bam_qc.sh $tissue

# sort file and generate bigWig files. 
bash -x ./sort_file.bigWig.sh $tissue

# count the reads into sparse matrix. 
# bash ./count2bin.sh


# call peaks. 
bash -x callPeaks.sh $tissue
# count the reads to peaks and construct matrix. 
bash -x count2peaks.sh $tissue

######## OUT OF DEPENDENCY #########

# calculate barcode 
Rscript plot_num_barcodes.r $tissue
