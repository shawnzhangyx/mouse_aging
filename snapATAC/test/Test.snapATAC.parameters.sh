##### Default parameter: ####

### Barcode selection. 
# TSS_enrich = 10 [0,7,10]
# fragments = 500 [500,1000]

### Dimension reduction.
# landmark sampleing: seed: 1 [1-5]
# landmark cells: 10000 [5000,10000,20000,40000]
# dim PC = 20  [5,10,20,50]

### Clustering: 
# k = 50 [10,30,50,100]
# res = 0.5 [0.1,0.3,0.5,0.8]
##########################


## Test barcode selection. 
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_barcodes/

for TSS in 0 7 10; do
  for frag in 500 1000; do
Rscript Test.snapATAC.cluster.TssEnrichment.R \
    -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
    --fragment-cutoff $frag \
    --tss-cutoff $TSS \
    -o ../../analysis/snapATAC/DH/snapFiles/Tests_barcodes/DH.pool.snapATAC &
  done
done


## Test dimension reduction. 
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_seeds/
for seed in {1..5}; do
Rscript Test.snapATAC.cluster.TssEnrichment.R \
    -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
    --seed $seed \
    -o ../../analysis/snapATAC/DH/snapFiles/Tests_seeds/DH.pool.snapATAC &
done
# landmark cells. 
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_landmarks/
for landmark in 5000 10000 20000 40000; do 
Rscript Test.snapATAC.cluster.TssEnrichment.R \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --landmark_num $landmark \
  -o ../../analysis/snapATAC/DH/snapFiles/Tests_landmarks/DH.pool.snapATAC
done
# PC numbers.
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_PC_num/
for pc in 5 10 20 50; do
Rscript Test.snapATAC.cluster.TssEnrichment.R \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --pc_dim $pc \
  -o ../../analysis/snapATAC/DH/snapFiles/Tests_PC_num/DH.pool.snapATAC 
  done



## Test clustering parameters. 
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_knn_k/
for k in 10 30 50 100; do 
Rscript Test.snapATAC.cluster.TssEnrichment.R \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --knn_k $k \
  -o ../../analysis/snapATAC/DH/snapFiles/Tests_knn_k/DH.pool.snapATAC
done
# res 
mkdir ../../analysis/snapATAC/DH/snapFiles/Tests_leiden_res/
for res in 0.1 0.3 0.5 0.8; do 
Rscript Test.snapATAC.cluster.TssEnrichment.R \
  -i ../../analysis/snapATAC/DH/snapFiles/DH.pool.snapATAC.raw.RData \
  --leiden_res $res \
  -o ../../analysis/snapATAC/DH/snapFiles/Tests_leiden_res/DH.pool.snapATA
done





