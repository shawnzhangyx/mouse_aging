import pysam

inBam = pysam.AlignmentFile(inFile)
outBam = pysam.AlignmentFile(barcode+".bam","wb",template=inBam)
for line in bamF.fetch(until_eof=True):
  

