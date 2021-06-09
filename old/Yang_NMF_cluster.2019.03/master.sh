# as a beginning, try different ranks. 
tissue=heart
RMIN=5
RMAX=20

bash -x try_different_ranks.sh $tissue $RMIN $RMAX


# generate bam and bigWig file for a chosen rank
rank=10
bash -x generate_bam_bigWig_rankR.sh $tissue $rank
### run each replicates separately
bash -x run_replicates_separate.sh $tissue $rank

# make plots 
bash -x generate_DimRed_plots.sh $tissue $rank
# output the cluster promoter and peak fold change. 
Rscript output_cluster_feature.edger.r $tissue $rank




# call cluster specific peaks. 
bash -x callPeaks_by_cluster.sh $tissue $rank 
Rscript cluster.difftest.by_cluster.r $tissue $rank

# look for differential peaks during aging. 
Rscript cluster.difftest.r $tissue $rank
# motif analysis of differential peaks during aging
bash -x motif_analysis_diff_peaks.sh $tissue $rank 

# great analysis of differential peaks during aging
bash -x great_analysis_diff_peaks.sh $tissue $rank 


