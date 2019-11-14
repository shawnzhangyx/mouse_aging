import sys
from scipy.sparse import load_npz
import numpy as np
import argparse 

parser = argparse.ArgumentParser(description='Run NMF using sklearn.')
parser.add_argument('-i', '--inputF', type=str, dest="inputF", help='input matrix in npz format')
parser.add_argument('-x', '--xgi', type=str, dest="xgi", help='input xgi index')
parser.add_argument('-o', '--outF', type=str, dest="outF", help='output file')

args = parser.parse_args()


V = load_npz(args.inputF)
with open(args.xgi, 'r') as f:
  barcodes = [x.strip() for x in f.readlines()]
counts = np.sum(V,1)

if len(barcodes) < len(counts):
  print("length of barcodes is shorter than length of counts. Exit!")
  exit(1)

#with open(args.outF,'w') as f:
for idx in range(len(counts)):
  print(barcodes[idx]+"\t"+str(int(counts[idx])))
