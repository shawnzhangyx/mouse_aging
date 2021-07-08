# FC specific
for age in 03 18; do 
  #L23
  samtools merge -f bam_merge_rep/L23.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C0.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C0.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C2.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C2.bam
  #L4 
samtools merge -f bam_merge_rep/L4.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C1.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C1.bam
  #L6 
  samtools merge -f bam_merge_rep/L6.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C5.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C5.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C7.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C7.bam
  #L5-Unknown
  samtools merge -f bam_merge_rep/L5Unknown.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C11.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C11.bam
  #L5-Parm1
  samtools merge -f bam_merge_rep/L5Parm1.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C14.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C14.bam
  #L5-Fezf2
  samtools merge -f bam_merge_rep/L5Fezf2.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C17.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C17.bam
  #L5a-Deptor
  samtools merge -f bam_merge_rep/L5Deptor.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C19.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C19.bam
  done

## HC specific 
for age in 03 18; do
  #DG
  samtools merge -f bam_merge_rep/DG.${age}.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C4.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C4.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C6.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C6.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C16.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C16.bam
  #CA1 
  samtools merge -f bam_merge_rep/CA1.${age}.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C9.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C9.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C12.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C12.bam
  #CA2/3
  samtools merge -f bam_merge_rep/CA23.${age}.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C15.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C15.bam
  #Claustram 
  samtools merge -f bam_merge_rep/Claustrum.${age}.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C20.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C20.bam
  #CP
  samtools merge -f bam_merge_rep/CP.${age}.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C22.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C22.bam
  done

## Both 
for age in 03 18; do 
  #InNeu-Pvalb
  samtools merge -f bam_merge_rep/InNeuPvalb.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C3.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C3.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C3.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C3.bam
  #ASC
  samtools merge -f bam_merge_rep/ASC.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C8.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C8.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C8.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C8.bam
  #OGC
  samtools merge -f bam_merge_rep/OGC.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C10.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C10.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C10.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C10.bam
  #InNeu-Sst
  samtools merge -f bam_merge_rep/InNeuSst.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C13.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C13.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C13.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C13.bam
  #InNeu-Vip
  samtools merge -f bam_merge_rep/InNeuVip.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C18.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C18.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C18.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C18.bam
  #OPC 
  samtools merge -f bam_merge_rep/OPC.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C21.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C21.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C21.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C21.bam
  #Endo
  samtools merge -f bam_merge_rep/Endo.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C23.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C23.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C23.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C23.bam
  #Mgc
  samtools merge -f bam_merge_rep/MGC.${age}.bam bam_split/split_M${age}_FC_rep1.bam_Ageing_H3K9me3/M${age}_FC_rep1.bam_C24.bam bam_split/split_M${age}_FC_rep2.bam_Ageing_H3K9me3/M${age}_FC_rep2.bam_C24.bam bam_split/split_M${age}_HC_rep1.bam_Ageing_H3K9me3/M${age}_HC_rep1.bam_C24.bam bam_split/split_M${age}_HC_rep2.bam_Ageing_H3K9me3/M${age}_HC_rep2.bam_C24.bam
  done

