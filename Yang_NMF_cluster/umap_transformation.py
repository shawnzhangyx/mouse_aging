import numpy as np
import umap
import pandas
#from plotnine import *
import sys

### Main
normH=sys.argv[1]
statH=sys.argv[2]
metric = sys.argv[3]
neighbors = int(sys.argv[4])
output = sys.argv[5]

normH = np.loadtxt(normH)
statH = np.genfromtxt(statH,dtype=None,names=True)

reducer = umap.UMAP(random_state=42, metric=metric,n_neighbors=neighbors, n_components=2)
embedding = reducer.fit_transform(normH.T)
df = pandas.DataFrame({"x":embedding[:,0],"y":embedding[:,1],"cluster":statH["class0"]})

np.savetxt(output, df)

