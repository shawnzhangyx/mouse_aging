from scipy.io import mmwrite
from scipy.sparse import load_npz
import sys
inFile = sys.argv[1]
outFile = sys.argv[2]
print(inFile,outFile)

# read the npz file
a = load_npz(inFile)
mmwrite(outFile,a)
