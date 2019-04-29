### environments and variables. 
tissue=$1
MIN_R=$2
MAX_R=$3

mat_path=../../data/snATAC/ct2peaks/
out_path=../../analysis/Yang_NMF_method/${tissue}/

mkdir $out_path

for rank in $(seq ${MIN_R} ${MAX_R}); do 
  echo $rank
  mkdir $out_path/R${rank}
python snATAC_YLee/snATAC.nmf.lite.py -i $mat_path/${tissue}.all.npz \
      -x $mat_path/${tissue}.all.xgi \
      -y $mat_path/${tissue}.all.ygi \
      -o $out_path/R${rank}/${tissue}.R${rank} \
      -r $rank -n 1 -p 0.05 -c 500 > $out_path/R${rank}/${tissue}.R${rank}.log &
done 
wait

for rank in $(seq $MIN_R $MAX_R); do
(
# plot H & W
#Rscript snATAC_YLee/snATAC.plotH.R -i $out_path/R$rank/${tissue}.R${rank}.H.mx -o $out_path/R$rank/${tissue}.R${rank}.H
#Rscript snATAC_YLee/snATAC.plotW.R -i $out_path/R$rank/${tissue}.R${rank}.W.mx -o $out_path/R$rank/${tissue}.R${rank}.W

python snATAC_YLee/snATAC.nmf.stat.py -m $out_path/R$rank/${tissue}.R${rank}.npz \
      -x $out_path/R$rank/${tissue}.R$rank.xgi \
      -y $out_path/R$rank/${tissue}.R$rank.ygi \
      --basis $out_path/R$rank/${tissue}.R${rank}.W.mx \
      --coef $out_path/R$rank/${tissue}.R${rank}.H.mx \
      -c 0.2 -o $out_path/R$rank/${tissue}.R${rank}

# calculate cell sparseness and entropy using the statH file
#Rscript snATAC_YLee/snATAC.statBox.R -i $out_path/${tissue}.R${rank}.statH \
#      -o $out_path/${tissue}.R${rank} >> $out_path/${tissue}.R${rank}.sta.txt

# calculate silhouette and plot tSNE using coefficient matrix H
#python snATAC_YLee/snATAC.nmf.plot.py --normH $out_path/${tissue}.R${rank}.normH \
#      --statH $out_path/R$rank/${tissue}.R${rank}.statH \
#      -p 20 -o $out_path/plots/${tissue}.R${rank}

) &
done
wait

# annotate age and cluster to cells. 
for rank in $(seq $MIN_R $MAX_R); do
Rscript annotate_age_and_cluster_to_cells.r $out_path/R$rank/$tissue.R$rank.statH $out_path/R$rank/$tissue.R$rank.barcode.cluster.stage.rep.txt
Rscript plot_celltype_fraction.r $out_path/R$rank/$tissue.R$rank.statH $out_path/R$rank/$tissue.R$rank.celltype_fraction.pdf

done


Rscript snATAC_YLee/combine_cluster_info.for_sankey_plot.r --name_col 1 --label_col 3 --output sankey_plot.html $(ls $out_path/R*/$tissue.R*.statH -tr) 
mv sankey_plot.html $out_path


