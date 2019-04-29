#!/usr/bin/env python
from scipy.sparse import load_npz,save_npz
import numpy as np
import sys

matF = sys.argv[1]
bcF = sys.argv[2]
outF = sys.argv[3]


mat = load_npz(matF)
rowSum = np.sum(mat,axis=1)

with open(bcF,'r') as f:
  barcodes = [x.strip() for x in f.readlines()]

fout = open(outF,'w')
for idx,barcode in enumerate(barcodes):
   fout.write(barcode+"\t"+str(rowSum[idx,0])+"\n")
fout.close()





