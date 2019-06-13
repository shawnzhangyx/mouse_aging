### obtain the names of the repeats. 
bash -x build_repeat_names.sh

### intersect the bam file with the repeat annotation. 
bash -x intersect_uniq_bam_w_rep.sh


### count the repeats to each cell. 
bash -x count2repeat.sh

