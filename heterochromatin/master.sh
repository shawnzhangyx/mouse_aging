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

#computeMatrix scale-regions -S ../../analysis/snapATAC/FC/bigWig.cluster_age/FC.*.03.*.bw \
#  -R L1MA5A.bed \
#  --beforeRegionStartLength 2000 \
#  --regionBodyLength 5000 \
#  --afterRegionStartLength 2000 \
#  --skipZeros -o cluster.03.matrix.mat.gz

#plotProfile -m cluster.03.matrix.mat.gz \
#              -out Profile.by.Cluster.03.pdf \
#              --numPlotsPerRow 1 \
#              --plotTitle "Test data profile"


