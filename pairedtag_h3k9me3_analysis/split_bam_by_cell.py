#!/bin/python

import argparse
from multiprocessing import Pool
import pysam


parser = argparse.ArgumentParser(description='filter bam based on QNAMES')
parser.add_argument('--file', type=str, dest="file", help='file')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')
args = parser.parse_args()


out_bam_dict = {}

bamF = pysam.AlignmentFile(args.file)

N = 0
for b in bamF.fetch(until_eof=True):
  N = N + 1
  if N % 10000 == 0:
    print(str(N) + " lines")
  barcode = "_".join(b.query_name.split(':')[7:11])
  if barcode in out_bam_dict:
        out_bam_dict[barcode].write(b)
  else:
    out_bam_dict[barcode] = pysam.AlignmentFile(args.outPrefix+barcode+".bam", "wb", template=bamF)
    out_bam_dict[barcode].write(b)


for key in out_bam_dict:
  out_bam_dict[key].close()






