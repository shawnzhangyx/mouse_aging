region="chr8\t20775078\t21969439"

tmpfile=$(mktemp)
echo -e  $region > $tmpfile

bedtools intersect -a ../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/peaks/dorsal_hippocampus.C3_peaks.narrowPeak -b $tmpfile > ../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/C3.defensin.peaks.bed

findMotifsGenome.pl ../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/C3.defensin.peaks.bed  mm10 ../../analysis/Yang_NMF_method/dorsal_hippocampus/R15/C3.defensin.peaks.homer
