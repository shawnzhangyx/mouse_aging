tissue=$1
rank=$2

out_path=../../analysis/Yang_NMF_method/${tissue}/
mkdir $out_path/R$rank/UMAP


#for neighbors in $(seq 15 15 200); do
#  echo $neighbors 
#  python test.umap_transformation.py euclidean $neighbors &
#  python test.umap_transformation.py correlation $neighbors &
#  done

neighbors=30
python test.umap_transformation.py $tissue $rank euclidean $neighbors &
python test.umap_transformation.py $tissue $rank correlation $neighbors &
wait
Rscript test.umap_plot.r $tissue $rank
