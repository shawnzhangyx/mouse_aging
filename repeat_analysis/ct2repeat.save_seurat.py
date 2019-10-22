from scipy import sparse
from scipy.io import mmwrite
import sys,os
import gzip 

reads = sys.argv[1]
xgi = sys.argv[3]
ygi = sys.argv[2]
X_CUT = int(sys.argv[4])
outPrev = sys.argv[5]

barcodes = {}
barcode_cnts = {}
peaks = {}

## make directory
try: os.mkdir(outPrev)
except FileExistsError:
  pass

## read and write barcodes as XGI. 
with open(xgi,'r') as in_xgi:
  for idx,line in enumerate(in_xgi):
    bc,count = line.strip().split()
    if int(count) > X_CUT: 
      barcodes[bc] = idx
      barcode_cnts[bc] = int(count)

with gzip.open(outPrev+"/barcodes.tsv.gz",'wb') as out_xgi:
  for key, val in barcodes.items():
    out_xgi.write((key + "\n").encode())

## read and write peaks as YGI
with open(ygi,'r') as in_ygi:
  for idx,line in enumerate(in_ygi):
      peaks["\t".join(line.strip().split("\t")[:3])] = idx
  peaks["Others"] = len(peaks)

with gzip.open(outPrev+"/features.tsv.gz",'w') as out_ygi:
  for key, val in peaks.items():
    out_ygi.write((key+"\t"+key+"\t"+"Gene Expression\n").encode())

## read and process the count2repeat entries. 
with gzip.open(reads, 'rt') as inFile:
  ct_dict = {}
  prev_qname = ""
  for idx,line in enumerate(inFile):
    if idx >1 and idx % 1000000 == 0:
      print("%s lines processed"%idx) 
#      break
#    print(prev_qname)
    items = line.split('\t')
    peak = items[3] #"\t".join(items[:3])
    if items[4].split("/")[0] == prev_qname:
      continue
    prev_qname = items[4].split("/")[0]
    barcode = items[4].split(":")[0]
    # skip a peak if the annotation is not present. 
    if peak not in peaks: 
      continue
    if barcode in barcodes: 
      barcode_cnts[barcode] -= 1
    #could add the counts up
      try: 
        ct_dict[(barcodes[barcode],peaks[peak])] += 1
      except KeyError: 
        ct_dict[(barcodes[barcode],peaks[peak])] = 1
#      xgi_list.append(barcodes[barcode])
#      ygi_list.append(peaks[peak])
  ct_list = list(ct_dict.values()) + list(barcode_cnts.values())
  xgi_list = [x[0] for x in ct_dict.keys()] + [barcodes[x] for x in list(barcode_cnts.keys())]
  ygi_list = [x[1] for x in ct_dict.keys()] + [len(peaks)-1] * len(barcode_cnts)
#  print(ct_list[:5])
#  print(xgi_list[:5])
#  print(ygi_list[:5])
  print(len(ct_list),len(xgi_list),len(ygi_list))
  cooMx = sparse.coo_matrix((ct_list, (ygi_list,xgi_list)))
  #print(cooMx)
  #  sparse.save_npz(outPrev+".npz", cooMx)
  mmwrite(outPrev+"/matrix",cooMx)
  os.system("gzip -f "+outPrev+"/matrix.mtx")
