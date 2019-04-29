tissue=dorsal_hippocampus

# run snap_pre.sh on each file separately
bash -x run_snaptools_pre.sh $tissue
# run snap_pre.merged.R to combine all sample in the same tissue. 
bash -x merge_and_cluster_one_tissue.sh $tissue

