### obtain the names of the repeats. 
bash -x build_repeat_names.sh


tissue=DH
### intersect the bam file with the repeat annotation. 
bash -x intersect_uniq_bam_w_rep.sh $tissue


### count the repeats to each cell. 
### currently using seurat pacakge, but could use other things as well. 
bash -x count2repeat.sh $tissue

### plot the classes of repeats during aging. 
Rscript analyze_repeat_seurat.allsamples.r $tissue
Rscript analyze_repeat_seurat.allsamples.post_analysis.r $tissue

