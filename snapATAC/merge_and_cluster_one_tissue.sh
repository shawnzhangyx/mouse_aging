tissue=$1
path=../../analysis/snapATAC/$tissue/
resolution=0.7

cd $path/snapFiles/
Rscript ../../../../scripts/snapATAC/snapATAC.pre.merged.R \
    -i $tissue.snap.list \
    --fragment_num 500 \
    --mito_ratio 0.2 \
    --dup_ratio 0.5 \
    --umap_ratio 0.8 \
    --pair_ratio 0.7 \
    --bin_size 5000 \
    --pc_num 50 \
    --black_list /projects/ps-renlab/mandyjiang/doublet/mm10.blacklist.v2.customized.bed \
    --cpu 5 \
    -o $tissue.pool.snapATAC

cd - 

Rscript snapATAC.cluster.R \
    -i $path/snapFiles/$tissue.pool.snapATAC.pre.RData \
    -d 20 \
    -o $path/snapFiles/$tissue.pool.snapATAC

Rscript annotate_age_and_cluster_to_cells.r \
  $path/snapFiles/$tissue.pool.snapATAC.cluster.meta.txt \
  $path/$tissue.pooled.barcode.cluster.stage.rep.txt 



