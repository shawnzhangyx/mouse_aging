tissue=$1
rank=$2

path=../../analysis/Yang_NMF_method/${tissue}/R$rank/age_diff_bycluster
mkdir -p $path/great/

##great analysis
Rscript great.make_job_list.r $path $tissue $rank
~/software/great/greatBatchQuery.py $path/great_jobs.txt

for file in $(ls $path/great/*great.tsv)
  do
  echo "cut -f 1-22 $file > ${file/.tsv/.red.tsv}"
  cut -f 1-22 $file > ${file/.tsv/.red.tsv}
  done


Rscript great.plot_great_results.r $tissue $rank


