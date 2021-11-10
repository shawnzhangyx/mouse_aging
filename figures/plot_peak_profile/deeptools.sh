file03=../../../analysis/snapATAC/*/bigWig.cluster_age/*.metacell_*.03.sorted.rpkm.bw
file10=../../../analysis/snapATAC/*/bigWig.cluster_age/*.metacell_*.10.sorted.rpkm.bw
file18=../../../analysis/snapATAC/*/bigWig.cluster_age/*.metacell_*.18.sorted.rpkm.bw

computeMatrix reference-point -S \
  $file03 \
  -R pisd_ps3.bed \
  --referencePoint center \
  --averageTypeBins median \
  --skipZeros -o m03.mat.gz

computeMatrix reference-point -S \
  $file10 \
  -R pisd_ps3.bed \
  --referencePoint center \
  --averageTypeBins median \
  --skipZeros -o m10.mat.gz

computeMatrix reference-point -S \
  $file18 \
  -R pisd_ps3.bed \
  --referencePoint center \
  --averageTypeBins median \
  --skipZeros -o m18.mat.gz


plotProfile -m m03.mat.gz \
  --perGroup \
  --plotType heatmap \
  -out m03.profile.pdf
#  --numPlotsPerRow 1 \

plotProfile -m m18.mat.gz \
  --perGroup \
  --plotType heatmap \
  -out m18.profile.pdf

