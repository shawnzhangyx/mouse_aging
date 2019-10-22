library("SnapATAC")
setwd("../../analysis/snapATAC/HT/")
load("snapFiles/HT.pool.snapATAC.cluster.RData")

rmDoublets <- function(obj, mat, ...) {
  UseMethod("rmDoublets", obj);
}

#' @export
rmDoublets.default <- function(
    obj, 
    mat = c("pmat", "bmat", "gmat"),
    path.to.python = NULL,
    path.to.condaenv = NULL,
    path.to.venv = NULL,
    expected.doublet.rate = 0.06,
    min.counts = 3, 
    min.cells = 3L, 
    min.gene_variability_pctl = 85, 
    n.prin_comps = 30L
){
  
  cat("Epoch: checking input parameters ... \n", file = stderr())
    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is.snap(obj)){
            stop("obj is not a snap object");
        }       
        if((x=nrow(obj))==0L){
            stop("obj is empty");
        }
    }
  
  ncell = nrow(obj);
        
    mat = match.arg(mat);
    mat.use = methods::slot(obj, mat);
    if((x=nrow(mat.use)) == 0L){
        stop("mat is empty, add matrix to snap object first")
    }
        
    if((x=nrow(mat.use)) != ncell){
        stop("mat has different length with cell number, re-add this matrix to snap object");
    }
  
  cat("Epoch: loading python environments and packages ... \n", file = stderr())
  # load library and python env
  library("reticulate")
  if(!is.null(path.to.python)){
    reticulate::use_python(path.to.python)
    cat("use the Python located in:", path.to.python, "\n")
  }else if(!is.null(path.to.condaenv)){
    reticulate::use_virtualenv(path.to.condaenv)
    cat("use the Python in Conda environment:", path.to.condaenv, "\n")
  }else if(!is.null(path.to.venv)){
    reticulate::use_virtualenv(path.to.venv)
    cat("use the Python in virtual environments:", path.to.venv, "\n")
  }else{
    cat("By default, uses the version of Python found on your PATH", "\n")
  }

  # convert paramter to python integer
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)
  mat.use = r_to_py(mat.use)
    min.cells = as.integer(min.cells);
    n.prin_comps = as.integer(n.prin_comps);
  
  # load python packages
  os <- import("os", convert = FALSE)
  scipy.io <- import("scipy.io", convert = FALSE)
  np <- import("numpy", convert = FALSE)
  scr <- import("scrublet", convert = FALSE)
  
  # identify potential doublets
  cat("Epoch: identify potential doublets ... \n", file = stderr())
  scrub = scr$Scrublet(mat.use, 
                       expected_doublet_rate=expected.doublet.rate)
  
  out = scrub$scrub_doublets(min_counts = min.counts, 
                             min_cells = min.cells, 
                             min_gene_variability_pctl = min.gene_variability_pctl, 
                             n_prin_comps = n.prin_comps)
  # convert to R explicitly at the end
  scrub <- py_to_r(scrub)
  out <- py_to_r(out)
  
  # summary the results
  cat("Epoch: summary doublets detection results ... \n", file = stderr())
  message("Automatically set threshold at doublet score = ", round(scrub$threshold_,2));
  message("Detected doublet rate = ", round(scrub$detected_doublet_rate_ * 100,2), "%");
  message("Estimated detectable doublet fraction = ", round(scrub$detectable_doublet_fraction_ * 100, 2), "%");
  message("Overall doublet rate:");
    message("\tExpected   = ", round(scrub$expected_doublet_rate * 100, 2), "%");
    message("\tEstimated  = ", round(scrub$overall_doublet_rate_ * 100, 2), "%");

  # output
  outList <- list("scrub" = scrub, "doublet_scores" = out[[1]], "predicted_doublets" = out[[2]]);
  
  return(outList);
  }

res <- rmDoublets(x.sp, mat="bmat", path.to.venv = "/mnt/silencer2/home/yangli/apps/anaconda3/envs/venv_scrublet/")


library("scales")
plotValue <- function(obj, values, viz.method, point.size, point.color, point.shape,  background.point, background.point.color, background.point.alpha, background.point.size, background.point.shape, low.value, high.value, down.sample, seed.use, plot.nrow, plot.ncol, pdf.file.name, pdf.height, pdf.width,...){
  UseMethod("plotValue", obj);
}

plotValue.default <- function(
    obj,
    values,
    viz.method=c("tsne", "umap"),
    point.size=0.5,
    point.color="red",
    point.shape=19,
    background.point=TRUE,
    background.point.color="grey",
    background.point.alpha=0.3,
    background.point.size=0.5,
    background.point.shape=19,
    low.value=0.0,
    high.value=1.0,
    down.sample=5000,
    seed.use=10,
    plot.nrow=3,
    plot.ncol=3,
    pdf.file.name=NULL,
    pdf.height=7,
    pdf.width=7,
    ...
){

    if(missing(obj)){
        stop("obj is missing")
    }else{
        if(!is(obj, "snap")){
            stop("obj is not a snap object")
        }
        ncell = nrow(obj);
        data.use = values;
        if((x=length(data.use)) == 0L){
            stop("value is empty, input a list of value")
        }
        if((x = length(data.use))!= ncell){
            stop("value has different number of rows from cell number, add value again")
        }
    }

    viz.method = match.arg(viz.method);
    viz.use = methods::slot(obj, viz.method);

    if((x = nrow(viz.use)) == 0L){
        stop("visulization matrix is empty, run runViz first")
    }

    if(down.sample < ncell){
        set.seed(seed.use);
        idx = sort(sample(seq(ncell), down.sample));
        data.use = data.use[idx];
        viz.use = viz.use[idx,];
    }

    
    if(!is.null(pdf.file.name)){
        if(file.exists(pdf.file.name)){
            warning("pdf.file already exists");
            file.remove(pdf.file.name);
        }else{
            if(!file.create(pdf.file.name)){
                stop("cannot create pdf.file, not a directory")
            }
            file.remove(pdf.file.name);
        }
        pdf(pdf.file.name,width=pdf.width,height=pdf.height);
    }

    y = data.use;
    if(background.point){
        plot(viz.use,
                  main="plotValue",
             col=scales::alpha(background.point.color, background.point.alpha),
             cex=background.point.size,
             pch=background.point.shape,
                yaxt='n',
                xaxt="n",
                xlab="",
                ylab="",
                ...
             );
        points(viz.use,
                    col=alpha(point.color, y),
                    cex=point.size,
                    pch=point.shape
                  );
    }else{
        plot(viz.use,
               main="plotValue",
                 col=alpha(point.color, y),
                 cex=point.size,
                 pch=point.shape,
                    yaxt='n',
                    xaxt="n",
                    xlab="",
                    ylab="",
                 ...
                 );
        }

    if(!is.null(pdf.file.name)){
        dev.off()
    }
    par(mfrow=c(1,1));
}

