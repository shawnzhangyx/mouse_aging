tissue=$1
rank=$2

cd ../../analysis/Yang_NMF_method/${tissue}/R$rank/age_diff_bycluster/
mkdir -p motif.homer

for file in *.bed;do 
  findMotifsGenome.pl $file mm10 motif.homer/${file/.bed/.homer} &
  done


