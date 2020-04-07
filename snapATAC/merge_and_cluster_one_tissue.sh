tissue=$1
path=../../analysis/snapATAC/$tissue/

cd $path/snapFiles/
Rscript ../../../../scripts/snapATAC/snapATAC.pre.createSnapObject.R \
    -i $tissue.snap.list \
    --tissue $tissue \
    --bin_size 5000 \
    --black_list /projects/ps-renlab/yanxiao/projects/mouse_aging/annotations/mm10-blacklist.v2.bed \
    --cpu 5 \
    -o $tissue.pool.snapATAC

cd - 

#Rscript snapATAC.cluster.R \
#    -i $path/snapFiles/$tissue.pool.snapATAC.raw.RData \
#    -d 20 \
#    --fragment-cutoff 500 \
#    --tss-cutoff 7 \
#    -o $path/snapFiles/$tissue.pool.snapATAC

## DH
Rscript snapATAC.cluster.v2.r \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --seed 1 \
  --tss-cutoff 10 \
  --marker_genes DH.marker_genes.txt \
  --fragment-cutoff 500 \
  --landmark_num 40000 \
  --pc_dim 20 \
  -o ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC
# All cells.
Rscript snapATAC.cluster.all_cells.r \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --seed 1 \
  --tss-cutoff 10 \
  --fragment-cutoff 500 \
  --pc_dim 20 \
  -o ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC
# Run harmony
Rscript snapATAC.cluster.harmony.r \
  ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.cluster.RData \
  ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7 

# FC

Rscript snapATAC.cluster.v2.r \
  -i ../../analysis/snapATAC/FC/snapFiles/FC.pool.snapATAC.raw.RData \
  --seed 1 \
  --tss-cutoff 10 \
  --marker_genes FC.marker_genes.txt \
  --fragment-cutoff 500 \
  --landmark_num 40000 \
  --pc_dim 20 \
  -o ../../analysis/snapATAC/FC/snapFiles/FC.pool.snapATAC
# RUn harmony
Rscript snapATAC.cluster.harmony.r \
  ../../analysis/snapATAC/FC/snapFiles/FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.cluster.RData \
  ../../analysis/snapATAC/FC/snapFiles/FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7

# HT 
Rscript snapATAC.cluster.v2.r \
  -i ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC.raw.RData \
    --seed 1 \
    --tss-cutoff 7 \
    --marker_genes HT.marker_genes.txt \
    --fragment-cutoff 500 \
    --landmark_num 20000 \
    --pc_dim 20 \
    -o ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC

Rscript snapATAC.cluster.all_cells.r \
  -i ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC.raw.RData \
  --seed 1 \
  --tss-cutoff 7 \
  --fragment-cutoff 500 \
  --pc_dim 20 \
  -o ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC
Rscript snapATAC.cluster.harmony.r \
  ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.cluster.RData \
  ../../analysis/snapATAC/HT/snapFiles/HT.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7

# LM 

Rscript snapATAC.cluster.v2.r \
  -i ../../analysis/snapATAC/LM/snapFiles/LM.pool.snapATAC.raw.RData \
    --seed 1 \
    --tss-cutoff 7 \
    --marker_genes HT.marker_genes.txt \
    --fragment-cutoff 500 \
    --landmark_num 40000 \
    --pc_dim 20 \
    -o ../../analysis/snapATAC/LM/snapFiles/LM.pool.snapATAC
Rscript snapATAC.cluster.harmony.r \
  ../../analysis/snapATAC/LM/snapFiles/LM.pool.snapATAC.Frag500.TSS7.Landmark40k.seed1.dimPC20.K20.res0.7.cluster.RData \
  ../../analysis/snapATAC/LM/snapFiles/LM.pool.snapATAC.Frag500.TSS7.Landmark40k.seed1.dimPC20.K20.res0.7


## make barcode information. 
Rscript make_barcode_information.r \
  DH \
  DH.pool.snapATAC.Frag500.TSS10.AllCells.seed1.dimPC20.K20.res0.7.harmony.meta.txt 

Rscript make_barcode_information.r \
  HT \
  HT.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.harmony.meta.txt

Rscript make_barcode_information.r \
  FC \
  FC.pool.snapATAC.Frag500.TSS10.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.meta.txt

Rscript make_barcode_information.r \
  LM \
  LM.pool.snapATAC.Frag500.TSS7.Landmark40k.seed1.dimPC20.K20.res0.7.harmony.meta.txt

