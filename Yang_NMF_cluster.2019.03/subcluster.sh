$tissue=dorsal_hippocampus
rank=15
cluster=3
MIN_R=2 #$2
MAX_R=5 #$3


mat_path=../../analysis/Yang_NMF_method/$tissue/R$rank/ct2peaks/
# create cluster. 
out_path=../../analysis/Yang_NMF_method/${tissue}/R$rank/subcluster/
mkdir $out_path

for rank2 in $(seq ${MIN_R} ${MAX_R}); do
  echo $rank2
  mkdir $out_path/R${rank2}
  python snATAC_YLee/snATAC.nmf.lite.py -i $mat_path/${tissue}.C$cluster.all.npz \
      -x $mat_path/${tissue}.C$cluster.all.xgi \
      -y $mat_path/${tissue}.C$cluster.all.ygi \
      -o $out_path/R${rank2}/${tissue}.R${rank2} \
      -r $rank2 -n 1 -p 0.05 -c 500 > $out_path/R${rank2}/${tissue}.C$cluster.R${rank2}.log  &
  done

wait

for rank2 in $(seq $MIN_R $MAX_R); do
(
# plot H & W
#Rscript snATAC_YLee/snATAC.plotH.R -i $out_path/R$rank2/${tissue}.C$cluster.R${rank2}.H.mx -o $out_path/R$rank2/${tissue}.R${rank2}.H
#Rscript snATAC_YLee/snATAC.plotW.R -i $out_path/R$rank2/${tissue}.C$cluster.R${rank2}.W.mx -o $out_path/R$rank2/${tissue}.R${rank2}.W

python snATAC_YLee/snATAC.nmf.stat.py -m $out_path/R$rank2/${tissue}.R${rank2}.npz \
      -x $out_path/R$rank2/${tissue}.R$rank2.xgi \
      -y $out_path/R$rank2/${tissue}.R$rank2.ygi \
      --basis $out_path/R$rank2/${tissue}.R${rank2}.W.mx \
      --coef $out_path/R$rank2/${tissue}.R${rank2}.H.mx \
      -c 0.2 -o $out_path/R$rank2/${tissue}.R${rank2}

) &
done
wait

for rank2 in $(seq $MIN_R $MAX_R); do
Rscript annotate_age_and_cluster_to_cells.r $out_path/R$rank2/$tissue.R$rank2.statH $out_path/R$rank2/$tissue.R$rank2.barcode.cluster.stage.rep.txt
Rscript plot_celltype_fraction.r $out_path/R$rank2/$tissue.R$rank2.statH $out_path/R$rank2/$tissue.R$rank2.celltype_fraction.pdf
done


Rscript snATAC_YLee/combine_cluster_info.for_sankey_plot.r --name_col 1 --label_col 3 --output sankey_plot.html $(ls $out_path/R*/$tissue.R*.statH -tr)
mv sankey_plot.html $out_path


