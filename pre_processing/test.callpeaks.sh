cd ../../data/snATAC/

rep1="bam.filter/DH_03_rep1/DH_03_rep1.filter.bam bam.filter/DH_10_rep1/DH_10_rep1.filter.bam bam.filter/DH_18_rep1/DH_18_rep1.filter.bam"
rep2="bam.filter/DH_03_rep2/DH_03_rep2.filter.bam bam.filter/DH_10_rep2/DH_10_rep2.filter.bam bam.filter/DH_18_rep2/DH_18_rep2.filter.bam"

source activate py27 && macs2 callpeak -t $rep1 $rep2 -n peaks/test.q0.1.DH. -g mm -f BAMPE --keep-dup all -q 0.1

source activate py27 && macs2 callpeak -t $rep1 $rep2 -n peaks/test.q0.2.DH. -g mm -f BAMPE --keep-dup all -q 0.2

macs2 callpeak -t $rep1 $rep2 -n peaks/test.p0.05.DH. -g mm -f BAMPE --keep-dup all -p 0.05

macs2 callpeak -t $rep1 $rep2 -n peaks/test.SE.p0.05.DH -g mm --keep-dup all -p 0.05


