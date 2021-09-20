mkdir ../../analysis/heterochromatin 


## plot L1MA5A profile. 
cut -f 6-8,10 ~/annotations/mm10/annotation/repeats/repName/L1MA5A.txt > ../../analysis/heterochromatin/L1MA5A.bed 

computeMatrix scale-regions -S ../../analysis/snapATAC/FC/bigWig.cluster_age/FC.*.bw \
  -R L1MA5A.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 2000 \
  --skipZeros -o cluster_age.matrix.mat.gz


plotProfile -m matrix.mat.gz \
              -out Profile.by.Cluster_Age.pdf \
              --numPlotsPerRow 3 \
              --plotTitle "Test data profile"

computeMatrix scale-regions -S ../../analysis/snapATAC/DH/bigWig.cluster_age/DH.*.bw \
  -R L1MA5A.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 2000 \
  --skipZeros -o DH.cluster_age.matrix.mat.gz

## plot IAPLTR3-int profile 

cut -f 6-8,10 ~/annotations/mm10/annotation/repeats/repName/IAPLTR3-int.txt > ../../analysis/heterochromatin/IAPLTR3-int.bed

computeMatrix scale-regions -S ../../analysis/snapATAC/FC/bigWig.cluster_age/FC.*.bw \
  -R IAPLTR3-int.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 2000 \
  --skipZeros -o FC.IAPLTR3-int.cluster_age.matrix.mat.gz

plotProfile -m FC.IAPLTR3-int.cluster_age.matrix.mat.gz \
              -out FC.IAPLTR3-int.Profile.by.Cluster_Age.pdf \
              --numPlotsPerRow 3 \
              --plotTitle "Test data profile"


## test IAPLTR3-int profile
cd ../../analysis/heterochromatin

awk -v OFS="\t" '{ print $1":"$2"-"$3,$1,$2,$3,"." }' IAPLTR3-int.bed > IAPLTR3-int.saf
awk -v OFS="\t" '{ print $1":"$2"-"$3,$1,$2,$3,"." }' L1MA5A.bed > L1MA5A.saf

files=$(ls ../../analysis/snapATAC/FC/bam.cluster_age/*.sorted.bam)
featureCounts -a IAPLTR3-int.saf -o FC.IAPLTR3-int.counts $files -F SAF -T 16
featureCounts -a L1MA5A.saf -o FC.L1MA5A.counts $files -F SAF -T 16


