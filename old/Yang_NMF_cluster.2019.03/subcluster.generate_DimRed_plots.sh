tissue=$1
rank1=$2
rank2=$3

out_path=../../analysis/Yang_NMF_method/${tissue}/R$rank1/subcluster/
mkdir $out_path/R$rank2/UMAP


#for neighbors in $(seq 15 15 200); do
#  echo $neighbors 
#  python test.umap_transformation.py euclidean $neighbors &
#  python test.umap_transformation.py correlation $neighbors &
#  done

neighbors=15
metric=euclidean
normH=$out_path/R$rank2/$tissue.R$rank2.normH
statH=$out_path/R$rank2/$tissue.R$rank2.statH
outFile=$out_path/R$rank2/UMAP/$tissue.R$rank2.$metric.$neighbors.txt
python umap_transformation.py $normH $statH $metric $neighbors $outFile
#wait
Rscript umap_plot.r $out_path/R$rank2/UMAP
