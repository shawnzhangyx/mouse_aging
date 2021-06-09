set -x 

## if the server is silencer, run interactively; if on TSCC, us job system. 
server=$(hostname)
rule=$1  ## first argument specific rules to run

if [ $server == "silencer.sdsc.edu" ]; then
  snakemake --snakefile Snakefile_snap_pre $rule -j 30 --ri
elif [ $server == "tscc-login1.sdsc.edu" ] || [ $server == "tscc-login2.sdsc.edu" ]; then
  snakemake $rule --snakefile Snakefile_snap_pre -p  -k -j 1000 --ri \
  --cluster "qsub -l nodes=1:ppn={threads} -N {rule} -q hotel -o pbslog/{rule}.{params.pbsName}.pbs.out -e pbslog/{rule}.{params.pbsName}.pbs.err" \
  --jobscript ../../scripts/pre_processing/jobscript.pbs --jobname "{rulename}.{jobid}.pbs" 2> >(tee -a snakemake.log >&2)


else
  echo -e "Invalide server option: $server"; exit 1;
fi


##### quick fix for the peaks. 

source activate py27
for tissue in DH FC HT LM; do 
for age in 03 10 18; do 
for rep in rep1 rep2; do
#(
#snapfile=../../analysis/snapATAC/$tissue/snapFiles/${tissue}_${age}_${rep}.snap
#echo $snapfile
#snaptools snap-del --session-name PM --snap-file=$snapfile
#echo add-pmat $snapfile 
#snaptools snap-add-pmat \
#    --snap-file=$snapfile  \
#    --peak-file=../../data/snATAC/peaks/${tissue}_summits.ext1k.filter_highSignal_repeats.bed
#) & 
done
done
done
