library(biomaRt)

setwd("../../annotations")
a=read.csv("genage_human.csv",stringsAsFactors=F)
b=read.csv("genage_models.csv",stringsAsFactors=F)
# Basic function to convert mouse to human gene names
convertHumanGeneList <- function(x){
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
mousex <- unique(genesV2[, 2])
return(mousex)
}

musGene = convertHumanGeneList(a$symbol)
musGene2 = c(musGene,b$symbol[which(b$organism == "Mus musculus")])

write.table(unique(musGene2),"genage.mouse.symbol.txt",col.names=F,quote=F,row.names=F)

