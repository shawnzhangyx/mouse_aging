# link files
cd ../../analysis/snapATAC/
mkdir shared.endo shared.mac
ln DH/age_diff_edgeR.snap/13.up.bed shared.endo/DH.up.bed
ln FC/age_diff_edgeR.snap/15.up.bed shared.endo/FC.up.bed
ln HT/age_diff_edgeR.snap/3.up.bed shared.endo/HT.3.up.bed
ln HT/age_diff_edgeR.snap/7.up.bed shared.endo/HT.7.up.bed
ln LM/age_diff_edgeR.snap/6.up.bed shared.endo/LM.up.bed
ln BM/age_diff_edgeR.snap/15.up.bed shared.endo/BM.up.bed

ln DH/age_diff_edgeR.snap/13.down.bed shared.endo/DH.down.bed
ln FC/age_diff_edgeR.snap/15.down.bed shared.endo/FC.down.bed
ln HT/age_diff_edgeR.snap/3.down.bed shared.endo/HT.3.down.bed
ln HT/age_diff_edgeR.snap/7.down.bed shared.endo/HT.7.down.bed
ln LM/age_diff_edgeR.snap/6.down.bed shared.endo/LM.down.bed
ln BM/age_diff_edgeR.snap/15.down.bed shared.endo/BM.down.bed

ln DH/age_diff_edgeR.snap/13.edger.txt shared.endo/DH.edger.txt
ln FC/age_diff_edgeR.snap/15.edger.txt shared.endo/FC.edger.txt
ln HT/age_diff_edgeR.snap/3.edger.txt shared.endo/HT.3.edger.txt
ln HT/age_diff_edgeR.snap/7.edger.txt shared.endo/HT.7.edger.txt
ln LM/age_diff_edgeR.snap/6.edger.txt shared.endo/LM.edger.txt
ln BM/age_diff_edgeR.snap/15.edger.txt shared.endo/BM.edger.txt


mkdir shared.endo/bigWig/
ln DH/bigWig.cluster_age/DH.metacell_13* shared.endo/bigWig/
ln FC/bigWig.cluster_age/FC.metacell_15* shared.endo/bigWig/
ln HT/bigWig.cluster_age/HT.metacell_3* shared.endo/bigWig/
ln HT/bigWig.cluster_age/HT.metacell_7* shared.endo/bigWig/
ln LM/bigWig.cluster_age/LM.metacell_6* shared.endo/bigWig/
ln BM/bigWig.cluster_age/BM.metacell_15* shared.endo/bigWig/

mkdir shared.endo/bam/
ln DH/bam.cluster_age_rep/DH.metacell_13.*.rep?.bam shared.endo/bam/
ln FC/bam.cluster_age_rep/FC.metacell_15.*.rep?.bam shared.endo/bam/
ln HT/bam.cluster_age_rep/HT.metacell_3.*.rep?.bam shared.endo/bam/
ln HT/bam.cluster_age_rep/HT.metacell_7.*.rep?.bam shared.endo/bam/
ln LM/bam.cluster_age_rep/LM.metacell_6.*.rep?.bam shared.endo/bam/
ln BM/bam.cluster_age_rep/BM.metacell_15.*.rep?.bam shared.endo/bam/

mkdir shared.endo/peaks/
ln DH/peak.cluster/DH.metacell_13.replicated_peaks.narrowPeak shared.endo/peaks
ln FC/peak.cluster/FC.metacell_15.replicated_peaks.narrowPeak shared.endo/peaks
ln HT/peak.cluster/HT.metacell_3.replicated_peaks.narrowPeak shared.endo/peaks
ln HT/peak.cluster/HT.metacell_7.replicated_peaks.narrowPeak shared.endo/peaks
ln LM/peak.cluster/LM.metacell_6.replicated_peaks.narrowPeak shared.endo/peaks
ln BM/peak.cluster/BM.metacell_15.replicated_peaks.narrowPeak shared.endo/peaks


#ln DH/age_diff_edgeR.snap/10.up.bed shared.mac/DH.up.bed
#ln FC/age_diff_edgeR.snap/9.up.bed shared.mac/FC.up.bed
#ln HT/age_diff_edgeR.snap/5.up.bed shared.mac/HT.up.bed
#ln LM/age_diff_edgeR.snap/9.up.bed shared.mac/LM.up.bed
#ln BM/age_diff_edgeR.snap/10.up.bed shared.mac/BM.up.bed

#ln DH/age_diff_edgeR.snap/10.down.bed shared.mac/DH.down.bed
#ln FC/age_diff_edgeR.snap/9.down.bed shared.mac/FC.down.bed
#ln HT/age_diff_edgeR.snap/5.down.bed shared.mac/HT.down.bed
#ln LM/age_diff_edgeR.snap/9.down.bed shared.mac/LM.down.bed
#ln BM/age_diff_edgeR.snap/10.down.bed shared.mac/BM.down.bed

#mkdir shared.mac/bigWig/
#ln DH/bigWig.cluster_age/DH.metacell_10* shared.mac/bigWig/
#ln FC/bigWig.cluster_age/FC.metacell_9* shared.mac/bigWig/
#ln HT/bigWig.cluster_age/HT.metacell_5* shared.mac/bigWig/
#ln LM/bigWig.cluster_age/LM.metacell_9* shared.mac/bigWig/
#ln BM/bigWig.cluster_age/BM.metacell_10* shared.mac/bigWig/

