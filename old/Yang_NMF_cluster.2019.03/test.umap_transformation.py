### loading packages
import numpy as np
import umap
import pandas
#from plotnine import *
import sys

### Main
tissue=sys.argv[1]
rank=sys.argv[2]
metric = sys.argv[3]
neighbors = int(sys.argv[4])


normH = np.loadtxt("../../analysis/Yang_NMF_method/"+tissue+"/R"+rank+"/"+tissue+".R"+rank+".normH")
statH = np.genfromtxt("../../analysis/Yang_NMF_method/"+tissue+"/R"+rank+"/"+tissue+".R"+rank+".statH",dtype=None,names=True)

reducer = umap.UMAP(random_state=42, metric=metric,n_neighbors=neighbors, n_components=2)
embedding = reducer.fit_transform(normH.T)
df = pandas.DataFrame({"x":embedding[:,0],"y":embedding[:,1],"cluster":statH["class0"]})

np.savetxt("../../analysis/Yang_NMF_method/"+tissue+"/R"+rank+"/UMAP/" + 
          tissue + ".R" + rank + "." + metric + "." +
            str(neighbors)+".txt", df)


