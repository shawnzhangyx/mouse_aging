#filter the most conservative set of peak clusters. 
Rscript peak_cluster.distribute.r
# reduce the set between cell types. 
> clustered_diff_peaks.merged.bed
bedtools merge -i <(grep Down clustered_diff_peak.bed) |awk -v OFS="\t" '{print $0,"Down"}' >> clustered_diff_peaks.merged.bed
bedtools merge -i <(grep Up clustered_diff_peak.bed) |awk -v OFS="\t" '{print $0,"Up"}' >> clustered_diff_peaks.merged.bed

bedtools cluster -i clustered_diff_peak.bed > clustered_diff_peaks.clustered.bed




# plot Manhattan plot of peak clusters with smooth gaussain values.
Rscript  test.plot_age_change_peak_location.r
intersectBed -a clustered_diff_peaks.clustered.bed -b FB.P0.H3K9me3.gt100k_domain.bed -c > clustered_diff_peaks.clustered.ovlpH3K9me3.bed
# plot the number of peak cluster overlap with heterochromatin. 
Rscript plot.num_peak_cluster.barplot.r
# plot the base of genome covered by H3K9me3 domain and percent ovlp with peaks.
## need to run the heterochromatin scripts before this.(./scripts/heterochromatin)
intersectBed -a FB.P0.H3K9me3.gt100k_domain.bed \
  -b clustered_diff_peaks.merged.bed -wao \
  > FB.P0.H3K9me3.gt100k_domain.ovlp_ATACpeaks.bed
#  -b ../../../analysis/heterochromatin/DH_FC.up.peaks.bed -c |\
#  intersectBed -a - -b ../../../analysis/heterochromatin/DH_FC.down.peaks.bed -c|\
#  intersectBed -a - -b clustered_diff_peaks.clustered.bed -c > FB.P0.H3K9me3.gt100k_domain.ovlp_ATACpeaks.bed

