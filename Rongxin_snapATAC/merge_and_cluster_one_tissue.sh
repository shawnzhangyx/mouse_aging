tissue=$1
path=../../analysis/Rongxin_snapATAC/$tissue/
resolution=0.7


Rscript snapATAC.pre.merged.R -i $path/snapFiles/snap.list \
    --fragment_num 500 \
    --mito_ratio 0.2 \
    --dup_ratio 0.5 \
    --umap_ratio 0.8 \
    --pair_ratio 0.7 \
    --bin_size 5000 \
    --pc_num 50 \
    --black_list mm10.blacklist.bed.gz \
    --cpu 5 \
    -o $path/snapFiles/$tissue.pool.snapATAC

Rscript snapATAC.cluster.R \
    -i $path/snapFiles/$tissue.pool.snapATAC.pre.RData \
    -r $resolution \
    -o $path/snapFiles/$tissue..pool.snapATAC

Rscript annotate_age_and_cluster_to_cells.r \
  $path/snapFiles/$tissue.pool.snapATAC.cluster.dims25.k15.resolution${resolution}.meta.txt \
  $path/$tissue.pooled.barcode.cluster.stage.rep.txt \


