
id=$1 #DH_03_rep1 
sample=$2 #JB_55

#cellranger count --id=JB_55 \
#  --transcriptome=/projects/ps-renlab/yanxiao/software/cellranger-3.0.2/refdata-cellranger-mm10-3.0.0 \
#  --fastqs=/projects/ps-renlab/yanxiao/projects/mouse_aging/data/snRNA/fastq/ \
#  --sample=JB_55 \
#  --project=JB_55 \
#  --expect-cells=5000 \
#  --localmem=50 \
#  --localcores 16 

mkdir -p ../../data/snRNA/cellranger.intron/ && cd ../../data/snRNA/cellranger.intron/


cellranger count --id=$id \
  --transcriptome=/projects/ps-renlab/yanxiao/software/cellranger-3.0.2/mm10.3.0.0_premrna \
  --fastqs=/projects/ps-renlab/yanxiao/projects/mouse_aging/data/snRNA/fastq/ \
  --sample=$sample \
  --project=$sample \
  --expect-cells=5000 \
  --localmem=32 \
  --localcores 8

