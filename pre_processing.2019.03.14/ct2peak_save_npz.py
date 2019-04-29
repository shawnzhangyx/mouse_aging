from scipy import sparse
from scipy.io import mmwrite
import sys
import gzip 

reads = sys.argv[1]
xgi = sys.argv[3]
ygi = sys.argv[2]
X_CUT = int(sys.argv[4])
outPrev = sys.argv[5]

barcodes = {}
peaks = {}

## read and write barcodes as XGI. 
with open(xgi,'r') as in_xgi:
  for idx,line in enumerate(in_xgi):
    bc,count = line.strip().split()
    if int(count) > X_CUT: 
      barcodes[bc] = idx

with open(outPrev+".xgi",'w') as out_xgi:
  for key, val in barcodes.items():
    out_xgi.write(key + "\n")

## read and write peaks as YGI
with open(ygi,'r') as in_ygi:
  for idx,line in enumerate(in_ygi):
      peaks["\t".join(line.strip().split("\t")[:3])] = idx

with open(outPrev+".ygi",'w') as out_ygi:
  for key, val in peaks.items():
    out_ygi.write(key+"\n")

## read and process the count2peak entries. 

with gzip.open(reads, 'rt') as inFile:
  xgi_list = []
  ygi_list = []
  ct_list = []
  prev_qname = ""
  for idx,line in enumerate(inFile):
    if idx % 1000000 == 0:
      print("%s lines processed"%idx) 
#      break
#    print(prev_qname)
    items = line.split('\t')
    peak = "\t".join(items[:3])
    if items[3].split("/")[0] == prev_qname:
      continue
    prev_qname = items[3].split("/")[0]
    barcode = items[3].split(":")[0]
    if barcode in barcodes: 
      ct_list.append(1)
      xgi_list.append(barcodes[barcode])
      ygi_list.append(peaks[peak])
    #except KeyError:
    #  pass
  #print(ct_list,xgi_list,ygi_list)
  cooMx = sparse.coo_matrix((ct_list, (xgi_list,ygi_list)))
  #print(cooMx)
  sparse.save_npz(outPrev+".npz", cooMx)

