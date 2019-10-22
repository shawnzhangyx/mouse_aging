### obtain the names of the repeats. 
cut -f 4 /projects/ps-renlab/yanxiao/software/RepEnrich2/refs/setup_folder_mm10/repnames.bed|sort -u > ../../analysis/repeat_analysis/repnames.txt
### build the repeat name, class and family. 
Rscript build_repeat_class_family.r

