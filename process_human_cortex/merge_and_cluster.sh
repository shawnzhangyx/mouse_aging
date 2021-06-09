tissue=FC
cd ../../analysis/human_cortex/FC/snapFiles/

Rscript ../../../../scripts/process_human_cortex/snapATAC.pre.createSnapObject.R \
    -i $tissue.snap.list \
    --tissue $tissue \
    --bin_size 5000 \
    --black_list /projects/ps-renlab/yanxiao/projects/mouse_aging/annotations/hg38-blacklist.v2.bed \
    --cpu 10 \
    -o $tissue.pool.snapATAC

cd -

Rscript snapATAC.cluster.all_cells.r \
  -i ../../analysis/human_cortex/FC/snapFiles/FC.pool.snapATAC.raw.RData \
  --seed 1 \
  --tss-cutoff 7 \
  --fragment-cutoff 500 \
  --pc_dim 20 \
  -o ../../analysis/human_cortex/FC/snapFiles/FC.pool.snapATAC

# harmonize
Rscript ../snapATAC/snapATAC.cluster.harmony.r \
  ../../analysis/human_cortex/FC/snapFiles/FC.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.cluster.RData \
    ../../analysis/human_cortex/FC/snapFiles/FC.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7

Rscript make_barcode_information.r FC \
  FC.pool.snapATAC.Frag500.TSS7.AllCells.seed1.dimPC20.K20.res0.7.harmony.meta.txt

