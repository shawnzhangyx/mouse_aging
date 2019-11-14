import pysam
import pandas as pd
from multiprocessing import Pool


# fragment cutoff. 0, 500, 1000
# TSS enrichment cutoof: 0,4,7,10

# F0-T0 F0-T4 F0-T7 F0-T10
# F500-T0 F500-T4 F500-T7 F500-T10
# F1000-T0 F1000-T4 F1000-T7 F1000-T10



## generate bam file based on barcode and filenames
def generate_bams(inputBamName, barcodes,outputBamName): 
  bamF = pysam.AlignmentFile(inputBamName)
  obam = pysam.AlignmentFile(outputBamName,'wb',template=bamF)
  for b in bamF.fetch(until_eof=True):
    if b.query_name.split(':')[0] in barcodes:
      obam.write(b)
  obam.close()
  bamF.close()
  print("Finished")

## generate the bam files according to the cutoffs
meta = pd.read_table("qc/bc_info_by_sample/DH_18_rep1.barcode.info.txt")


F_CUT = 1000
TSS_CUT = 7
p = Pool(12)

for F_CUT in [0,500,1000]:
  for TSS_CUT in [0,4,7,10]:
    barcodes = meta[(meta["filter"] >= F_CUT) & (meta.TSS_enrich >=TSS_CUT)]["barcodes"]
    p.apply_async(generate_bams,("bam.filter.sort/DH_18_rep1/DH_18_rep1.filter.csort.bam",list(barcodes),"test.bam.cutoff/DH_18_rep1.F"+str(F_CUT)+"_TSS"+str(TSS_CUT)+".bam"))
    
p.close()
p.join()

