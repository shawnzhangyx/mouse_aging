> all.diff.bed
#for tissue in DH FC HT LM; do 
for tissue in DH FC ; do
  files=$(ls ../../../analysis/snapATAC/$tissue/age_diff_edgeR.snap/*.bed)
  cat $files|awk -v OFS="\t" -v t=$tissue '{print $1,$2,$3,t}'
  done | bedtools sort > all.diff.bed
bedtools merge -i all.diff.bed > merged.bed
intersectBed -a  merged.bed -b all.diff.bed -c > overlaped.c.bed
intersectBed -a  merged.bed -b all.diff.bed -wo > overlaped.wo.bed

#Rscript plot_num_overlaps.DH_FC.r

> all.diff.bed
for tissue in DH FC HT LM BM; do
  files=$(ls ../../../analysis/snapATAC/$tissue/age_diff_edgeR.snap/*.{up,down}.bed)
  cat $files|awk -v OFS="\t" -v t=$tissue '{print $1,$2,$3,t}'
  done | bedtools sort > all.diff.bed
bedtools merge -i all.diff.bed > merged.bed
intersectBed -a  merged.bed -b all.diff.bed -c > overlaped.c.bed
intersectBed -a  merged.bed -b all.diff.bed -wo > overlaped.wo.bed

Rscript plot_num_overlaps.r


### calculate the overlap by class
Rscript combine_class.r 
for class in ExN InN Glia Lymphoid Muscle Myeloid Other; do 
bedtools sort -i class_diff_peak/$class.bed > class_diff_peak/$class.sorted.bed
bedtools merge -i class_diff_peak/$class.sorted.bed |awk -v OFS="\t" -v c=$class '{print $0,c}' > class_diff_peak/$class.merged.bed 
done

cat class_diff_peak/*.merged.bed  |bedtools sort > class_diff_peak/All.sorted.bed
intersectBed -a  merged.bed -b class_diff_peak/ExN.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/Glia.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/InN.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/Lymphoid.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/Muscle.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/Myeloid.merged.bed -c |\
  intersectBed -a - -b class_diff_peak/Other.merged.bed -c \
  -c > overlaped.clades.bed

Rscript plot_num_overlaps.class.r
