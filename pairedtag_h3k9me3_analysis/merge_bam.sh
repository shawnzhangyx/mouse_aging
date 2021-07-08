# FC specific
for age in 03 18; do 
  for rep in rep1 rep2; do 
  #L23
  samtools merge -f bam_merge/L23.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C0.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C2.bam 
  #L6
  samtools merge -f bam_merge/L6.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C5.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C7.bam

  #L4 
  ln -s ../bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C1.bam bam_merge/L4.${age}.${rep}.bam
  #L5-Unknown
  ln -s ../bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C11.bam bam_merge/L5Unknown.${age}.${rep}.bam
  #L5-Parm1
  ln -s ../bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C14.bam bam_merge/L5Parm1.${age}.${rep}.bam 
  #L5-Fezf2
  ln -s ../bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C17.bam bam_merge/L5Fezf2.${age}.${rep}.bam 
  #L5a-Deptor
  ln -s ../bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C19.bam bam_merge/L5Deptor.${age}.${rep}.bam 
  done
done

## HC specific 
for age in 03 18; do
  for rep in rep1 rep2; do 
  #DG
  samtools merge -f bam_merge/DG.${age}.${rep}.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C4.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C6.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C16.bam
  #CA1 
  samtools merge -f bam_merge/CA1.${age}.${rep}.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C9.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C12.bam 
  #CA2/3
  ln -s ../bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C15.bam bam_merge/CA23.${age}.${rep}.bam
  #Claustrum 
  ln -s ../bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C20.bam bam_merge/Claustrum.${age}.${rep}.bam
  #CP
  ln -s ../bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C22.bam bam_merge/CP.${age}.${rep}.bam
  done
done

## Both 
for age in 03 18; do 
  for rep in rep1 rep2;do 
  #InNeu-Pvalb
  samtools merge -f bam_merge/InNeuPvalb.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C3.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C3.bam 
  #ASC
  samtools merge -f bam_merge/ASC.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C8.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C8.bam 
  #OGC
  samtools merge -f bam_merge/OGC.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C10.bam  bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C10.bam 
  #InNeu-Sst
  samtools merge -f bam_merge/InNeuSst.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C13.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C13.bam 
  #InNeu-Vip
  samtools merge -f bam_merge/InNeuVip.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C18.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C18.bam 
  #OPC 
  samtools merge -f bam_merge/OPC.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C21.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C21.bam 
  #Endo
  samtools merge -f bam_merge/Endo.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C23.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C23.bam 
  #Mgc
  samtools merge -f bam_merge/MGC.${age}.${rep}.bam bam_split/split_M${age}_FC_${rep}.bam_Ageing_H3K9me3/M${age}_FC_${rep}.bam_C24.bam bam_split/split_M${age}_HC_${rep}.bam_Ageing_H3K9me3/M${age}_HC_${rep}.bam_C24.bam 
  done
done
