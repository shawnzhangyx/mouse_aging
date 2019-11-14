library(dplyr)
library(Seurat)

pbmc.data <- Read10X(data.dir = "../../../analysis/repeat_analysis/DH/DH_03_rep1")
pbmc1 <- CreateSeuratObject(counts = pbmc.data, project = "DH", min.cells = 1, min.features = 1)
#pbmc <- NormalizeData(object = pbmc)

rsum1 = Matrix::rowSums(pbmc1)
rsum1[grep("rRNA",names(rsum1))]

#rep1
LSU-rRNA-Hsa SSU-rRNA-Hsa
       41305        26756
rsum1["GSAT-MM"]
 570376

> sum(rsum)
[1] 38360287
sum(a$UM[a$sample=="DH_03_rep1"],na.rm=T)
[1] 22340154


pbmc.data <- Read10X(data.dir = "../../../analysis/repeat_analysis/DH/DH_03_rep2")
pbmc2 <- CreateSeuratObject(counts = pbmc.data, project = "DH", min.cells = 1, min.features = 1)
#pbmc <- NormalizeData(object = pbmc)

rsum2 = Matrix::rowSums(pbmc2)

rsum2[grep("rRNA",names(rsum2))]


#rep2 
LSU-rRNA-Hsa SSU-rRNA-Hsa
       80152        55102
rsum2["GSAT-MM"]
GSAT-MM
 915914

sum(rsum)
[1] 115277002
sum(a$UM[a$sample=="DH_03_rep2"],na.rm=T)
[1] 43167599

