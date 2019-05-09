tissue=frontal_cortex
tissue=dorsal_hippocampus
tissue=heart

mat_path=../../data/snATAC/ct2peaks/
out_path=../../analysis/Yang_NMF_method/${tissue}/
mkdir $out_path

for rep in rep1 rep2; do 
(
mkdir $out_path/$rep

python snATAC_YLee/snATAC.nmf.lite.py -i $mat_path/${tissue}.$rep.npz \
      -x $mat_path/${tissue}.$rep.xgi \
      -y $mat_path/${tissue}.all.ygi \
      -o $out_path/$rep/${tissue} \
      -r 5 -n 1 -p 0.05 -c 500 > $out_path/$rep/${tissue}.log

# plot H & W.
#Rscript snATAC_YLee/snATAC.plotH.R -i $out_path/${tissue}.H.mx -o $out_path/${tissue}.H
#Rscript snATAC_YLee/snATAC.plotW.R -i $out_path/${tissue}.W.mx -o $out_path/${tissue}.W

python snATAC_YLee/snATAC.nmf.stat.py -m $out_path/${tissue}.npz \
      -x $out_path/$rep/${tissue}.xgi \
      -y $out_path/$rep/${tissue}.ygi \
      --basis $out_path/$rep/${tissue}.W.mx \
      --coef $out_path/$rep/${tissue}.H.mx \
      -c 0.2 -o $out_path/$rep/${tissue}

# calculate cell sparseness and entropy using the statH file
#Rscript snATAC_YLee/snATAC.statBox.R -i $out_path/${tissue}.statH \
#      -o $out_path/${tissue} >> $out_path/${tissue}.sta.txt

# calculate silhouette and plot tSNE using coefficient matrix H
python snATAC_YLee/snATAC.nmf.plot.py --normH $out_path/$rep/${tissue}.normH \
      --statH $out_path/$rep/${tissue}.statH \
      -p 20 -o $out_path/$rep/${tissue}
) &
done

