#!/usr/bin/env Rscript

sessionInfo()
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--tissue",required=TRUE, help="Tissue name (DH,FC,HT,LM)")
parser$add_argument("-i", "--input", required=TRUE, help="input txt with the information of path of snap file and name")
parser$add_argument("--bin_size", default = 5000, help="binSize to use [default %(default)s]")
parser$add_argument("--black_list", required=TRUE, help="black list file")
parser$add_argument("--cpu", default = 4, help="# of cpus [default %(default)s]")
parser$add_argument("-o", "--output", required=TRUE, help="output file prefix")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

suppressPackageStartupMessages(library("SnapATAC"))
library(tictoc)

inputF = args$input
bin_size = as.numeric(args$bin_size)
black_list = "/projects/ps-renlab/yanxiao/projects/mouse_aging/annotations/hg38-blacklist.v2.bed" #args$black_list
cpus =  as.numeric(args$cpu)
tissue= args$tissue
outF = args$output

inputf <- read.table(inputF,sep="\t",header=F)
fnames <- as.character(inputf[,2])
pnames <- as.character(inputf[,1])

tic("createSnap")
x.sp = createSnap(
    file=fnames,
    sample=pnames,
    num.cores=cpus
    );

x.sp = addBmatToSnap(
    x.sp,
    bin.size=bin_size,
    num.cores=cpus
    );
toc()

pdf(paste(outF, ".raw.sta.pdf",sep=""))
plotBarcode(x.sp, 
              pdf.file.name=NULL, 
              pdf.width=7, 
              pdf.height=7, 
              col="grey",
              border="grey",
              breaks=50
              );
dev.off()

pdf(paste(outF, ".pre.sta.pdf",sep=""))
plotBarcode(x.sp,
              pdf.file.name=NULL,
              pdf.width=7,
              pdf.height=7,
              col="grey",
              border="grey",
              breaks=50
              );
dev.off()

summarySnap(x.sp)

## Matrix Binarization
tic("makeBinary")
x.sp = makeBinary(x.sp, mat="bmat");
toc()

tic("calBmatCor")
calBmatCor(x.sp);
toc()

## Feature Selection (filter bins according to black list and chrom)
tic("filterBins")
black_list = read.table(black_list,sep="\t");
suppressPackageStartupMessages(library(GenomicRanges));
black_list.gr = GRanges(
                          black_list[,1],
                          IRanges(black_list[,2], black_list[,3])
                         );
idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr));
idy2 = grep("chrM|random", x.sp@feature);
idy = unique(c(idy1, idy2));
x.sp = x.sp[,-idy, mat="bmat"];

pdf(paste(outF, ".pre.BinCoverage.pdf",sep=""))
plotBinCoverage(
    x.sp,
    pdf.file.name=NULL,
    col="grey",
    border="grey",
    breaks=10,
    xlim=c(-6,6)
    );
dev.off()

x.sp = filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="bmat");
# filter empty rows
idx <- which(Matrix::rowSums(x.sp@bmat)!=0)
x.sp@metaData <- x.sp@metaData[idx,]
x.sp@barcode <- x.sp@barcode[idx]
x.sp@sample <- x.sp@sample[idx]
x.sp@file <- x.sp@file[idx]
x.sp@bmat <- x.sp@bmat[idx,]
toc()

## add the TSS enrichment information: 
x.sp@metaData$sample = x.sp@sample
meta_Files = list.files(pattern=tissue,
    path=paste0("../../../../data/human_cortex/qc/bc_info_by_sample"),
    full.names=T)

dat = list()
for (file in meta_Files){
  sample = sub(".*\\/(.._..._rep.).barcode.info.txt","\\1",file)
  dat[[sample]] =read.delim(file)
  dat[[sample]]$sample = sample
  }
dat2 = do.call(rbind,dat)
x.sp@metaData$TSS_enrich = dat2$TSS_enrich[match(paste(x.sp@metaData$sample,as.character(x.sp@metaData$barcode)),paste(dat2$sample, dat2$barcodes))]


## add peak matrix to the object: 
x.sp = addPmatToSnap(x.sp)


outfname = paste(outF, ".raw.RData",sep="")
save(x.sp, file=outfname)

outmeta = paste(outF, ".raw.meta.txt",sep="")
write.table(x.sp@metaData, file=outmeta, sep="\t", quote=F, col.names=T, row.names=F)


