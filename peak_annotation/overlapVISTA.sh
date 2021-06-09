
mkdir -p ../../analysis/peak_annotation/VISTA.ovlp/

for tissue in DH FC HT LM; do 
intersectBed -a /mnt/silencer2/home/shz254/datasets/vista_database/VISTA.hs_mm_combined.mm10.bed -b ../../data/snATAC/peaks/${tissue}_summits.ext1k.filter_highSignal_repeats.bed -wo > ../../analysis/peak_annotation/VISTA.ovlp/${tissue}.ovlp.VISTA.txt
done

mkdir -p ../../analysis/peak_annotation/VISTA.ovlp.diffPeak
for tissue in DH FC HT LM; do
intersectBed -a /mnt/silencer2/home/shz254/datasets/vista_database/VISTA.hs_mm_combined.mm10.bed -b <(cat ../../analysis/snapATAC/$tissue/age_diff_edgeR.snap/{?,??}.{up,down}.bed) -wo > ../../analysis/peak_annotation/VISTA.ovlp.diffPeak/${tissue}.ovlp.VISTA.txt
done
