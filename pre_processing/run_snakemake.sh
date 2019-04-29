set -x 

## if the server is silencer, run interactively; if on TSCC, us job system. 
server=$(hostname)
rule=$1  ## first argument specific rules to run

if [ $server == "silencer.sdsc.edu" ]; then
  snakemake $rule -j 30 
elif [ $server == "tscc-login1.sdsc.edu" ] || [ $server == "tscc-login2.sdsc.edu" ]; then
  snakemake $rule -p  -k -j 1000 --ri \
  --cluster "qsub -l nodes=1:ppn={threads} -N {rule} -q hotel -o pbslog/{rule}.{params.pbsName}.pbs.out -e pbslog/{rule}.{params.pbsName}.pbs.err" \
  --jobscript ../../scripts/pre_processing/jobscript.pbs --jobname "{rulename}.{jobid}.pbs" 

else
  echo -e "Invalide server option: $server"; exit 1;
fi
