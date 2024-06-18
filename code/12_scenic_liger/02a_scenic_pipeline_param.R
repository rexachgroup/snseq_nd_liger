#!/usr/bin/env Rscript
# Andrew Elkins
# 2019-07-16
# scenic_pipeline.R
# Running the scenic pipeline on our hoffman2 cluster

### Notes
# min genes and min samples set the default
# Need python for running GRNboost
# Need to have cisTarget database (feather) downloaded and use cis_path to point to directory
# For 40k cells, I used 128G of memory over 16 cores

### For Hoffman2
# Make sure to load:
# module load R/3.5.0
# module python/3.7.2
# Rscript scenic_pipeline.R expression_matrix_rds cell_info_rds --cis_path --nCores 10 --title ScenicSCBrain 

## 2020-02-10
#   Updates:
#     • Add ability to include own gene sets rather than genesKept
#     • Add argument with tsv file of genes to use

### Requirements -----------------------------------------------------------------------------
library(argparse)
suppressPackageStartupMessages(library(AUCell))
# suppressPackageStartupMessages(library(Seurat))
library(doParallel)
library(SCENIC)
library(reticulate)
# May or may not need this for running on hoffman, 
#reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)
reticulate::use_condaenv("../../conda/")

options(stringsAsFactors=FALSE)

### Functions --------------------------------------------------------------------------------
getArgs <- function(){
  parser <- ArgumentParser(description='Rscript to run SCENIC on human single cell data')
  # raw count matrix saved as an rds file
  parser$add_argument('expression_mat_rds', help=paste0('Path to expression matrix rds file (genes by cells)',
                                  ' Genes as rownames and Cells as colnames'),
                      default=NULL)
  # cell info has cells as rownames and a cluster number/name column
  parser$add_argument('cell_info_rds', 
                      help=paste0('Path to cell info rds file (clusters with cells as rownames)'),
                      default=NULL)
  parser$add_argument('out_dir', help=paste0("Folder to cd to. SCENIC int/ and output/ folders will be placed here."))
  parser$add_argument('--cis_path', 
                      help=paste0('Path to cis Target pathway directory(e.g. ~/home/scenic/cisTarget_databases)'),
                      default='cisTarget_databases')
  parser$add_argument('--nCores', 
                      help=paste0('Number of cores to use'), type="integer", default=1)
  parser$add_argument('--title', help=paste0('add project title name'),default="SCENIC Exp")
  
  parser$add_argument('--min_counts_gene', help=paste0('Filter by the total number of reads per gene. ',
                                                       'By default it keeps only the genes with at least',
                                                       ' 6 UMI counts across all samples '),
                      default=NULL)
  parser$add_argument('--min_samples', 
                      help=paste0('Filter by the number of cells in which the gene is detected'),
                      default=NULL)
  
  parser$add_argument('--genes_file',
                      help=paste0('The set of genes you would like SCENIC to use in the analysis instead of built in filtering. Should be tsv file.'),
                      default = NULL
                      )
  parser$parse_args()
}

createDirs <- function(out_dir){
    dir.create(file.path(out_dir, "output"), recursive = T)
    dir.create(file.path(out_dir, "int"), recursive = T)
}

### Procedure --------------------------------------------------------------------------------
main <- function(args){
  # Run through main routine
  if(is.null(args)){
    stop('Need to pass arguments!')
  } 
  if (!dir.exists(args$out_dir)) {
    dir.create(args$out_dir)
  }
  out_dir <- normalizePath(args$out_dir)
  cis_path <- normalizePath(args$cis_path)
  expression_mat_rds <- normalizePath(args$expression_mat_rds)
  cell_info_rds <- normalizePath(args$cell_info_rds)
  genes_file <- normalizePath(args$genes_file)
  
  get.pyfile <- normalizePath(list.files(pattern = "grn_boost.py$", recursive = TRUE, full.names = TRUE))
  print(get.pyfile)
  if(is.null(get.pyfile)){
    stop("Cant find grn_boost.py!")
  }

  print(out_dir)
  createDirs(out_dir)
  message(paste0(" chdir ", out_dir))
  setwd(out_dir)
  
  message("\nImport references...")
  # Fix for improperly-named motif annotation -- https://github.com/aertslab/SCENIC/issues/364#issuecomment-1366833415
  data(defaultDbNames)
  try({
    data(list = "motifAnnotations_hgnc_v9", package = "RcisTarget")
    .GlobalEnv$motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
  })
  dbs <- defaultDbNames[["hgnc"]]
  
  s.title <- ifelse(!is.null(args$title), "Scenic Analysis", args$title)

  scenic.options <- initializeScenic(
    org="hgnc", 
    dbDir=normalizePath(cis_path),
    dbs=dbs,
    datasetTitle=s.title, 
    nCores=args$nCores
  )
  message("\nDone import references")
  
  message("\nInitializing...")
  message("  ...reading in files")
  expr.mat <- readRDS(file=expression_mat_rds)
  cell.info <- readRDS(file=cell_info_rds)
  # Catch and formating
  if(ncol(expr.mat) != nrow(cell.info)){
    stop('Number of cells in two files does not match!')
  }
  # To fix cell names, sometimes have a "." instead of a "-"
  if(grepl("\\.",colnames(expr.mat)[1])){
    colnames(expr.mat) <- sapply(colnames(expr.mat), function(x) gsub("\\.", "-", x))
  }
  if(grepl("\\.",rownames(cell.info)[1] )){
    rownames(cell.info) <- sapply(rownames(cell.info), function(x) gsub("\\.", "-", x) )
  }
  if( !all.equal(colnames(expr.mat),rownames(cell.info)) ){
    stop('cell names from expression matrix and cell info do not match!')
  }


  colnames(cell.info)[1] <- "cluster"
  cell.info$nGene <- apply(expr.mat, 2, function(c) sum(c!=0))
  cell.info$nUMI <- apply(expr.mat, 2, sum)
  saveRDS(cell.info, file="int/cellInfo.Rds") # save cell info for later

  
  scenic.options@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  saveRDS(scenic.options, file="int/scenicOptions.Rds")

  message("  Saving scenic options")
  # Gene filtering step
  if (!is.null(genes_file)) {
    # If using own filtering on genes this should be true
    message("\tUsing given set of genes")
    genes.kept <- as.character(read.table(file = genes_file, sep = "\t")[,1])
    message(length(genes.kept), " genes")
  } else {
  # If min_counts_gene and min_samples not set, default to the SCENIC recomended values
    if (!is.null(args$min_counts_gene)) {
      min.counts.gene <- args$min_counts_gene
    } else {
      min.counts.gene <- 3*0.1*ncol(expr.mat)
    }
    if (!is.null(args$min_samples)) {
      min.samples <- args$min_samples
    } else {
      min.samples <- ncol(expr.mat) * 0.1
    }
    message("\tMin counts gene: ", min.counts.gene)
    message("\tMin samples: ", min.samples)
    
    message("\n...Filtering genes")
    genes.kept <- geneFiltering(
      expr.mat, 
      scenicOptions=scenic.options,
      minCountsPerGene=min.counts.gene,
      minSamples=min.samples
    )
    
    message("\tTotal kept genes: ", length(genes.kept))

  }
  
  saveRDS(genes.kept, file="int/1.1_genesKept.Rds")
  # Creating correlation matrix
  message("\n   Running correlation matrix")
  expr.mat.filtered <- expr.mat[rownames(expr.mat) %in% genes.kept, ]
  corr.mat <- cor(t(expr.mat.filtered), method="spearman")
  saveRDS(corr.mat, file=getIntName(scenic.options, "corrMat"))

  
  # Prep for GRNBoost
  message("   exporting for GRNBoost!")
  exportsForArboreto(expr.mat.filtered, scenicOptions=scenic.options)
  message("\nGRNBoost")
  run_arboreto(out_dir, get.pyfile) 
  #system("python", args = c(get.pyfile, out_dir))
  
  ## A fix since not using GENIE3, 
  # maybe make this modular 
  ####################################################################################
  # Update the file name and column name of GRNBoost output to match GENIE3 output
  grnboost.output <- importArboreto("int/1.2_grnoutput.txt")
  grnboost.output <- grnboost.output[, c("TF", "Target", "Importance")]
  colnames(grnboost.output) <- c("TF", "Target", "weight")
  saveRDS(grnboost.output, file="int/1.4_GENIE3_linkList.Rds")
  #####################################################################################
  
  # Standard SCENIC procedure, see link
  # https://rawcdn.githack.com/aertslab/SCENIC/0a4c96ed8d930edd8868f07428090f9dae264705/inst/doc/SCENIC_Running.html
  message("\nGenerating coexpression newtork modules")
  runSCENIC_1_coexNetwork2modules(scenic.options)
  
  message("\nGenerating regulons")
  runSCENIC_2_createRegulons(scenic.options)
  
  message("\nScoring cells")
  expr.mat.log <- log2(expr.mat + 1)
  rm(expr.mat)
  runSCENIC_3_scoreCells(
    scenic.options, 
    as.matrix(expr.mat.log), 
    skipTsne=TRUE,
    skipHeatmap=TRUE, 
    skipBinaryThresholds=TRUE
  )
  
  # Make some space
  rm(expr.mat.log)
  gc(verbose=FALSE)
  # Run binarization matrix procedure manually, for some reason the scenic 
  # wrapper was failing
  # make default thresholds (binarize matrix)
  cells.auc.thres <- NULL
  regulon.auc <- readRDS(file=getIntName(scenic.options, "aucell_regulonAUC"))
  cells.auc.thres <- AUCell_exploreThresholds(
    regulon.auc,
    smallestPopPercent=getSettings(scenic.options,"aucell/smallestPopPercent"),
    assignCells=TRUE, plotHist=FALSE,
    verbose=FALSE, nCores=scenic.options@settings$nCores
  )
  
  saveRDS(cells.auc.thres, file=getIntName(scenic.options, "aucell_thresholds"))
  message("\nProcess Complete!")
  
}

run_arboreto <- function(out_dir, grn_file) {
  # can't pass params to grn_boost.py through py_run_file -- set RSCENIC_PATH output dir
  Sys.setenv("RSCENIC_PATH" = out_dir)
  py_run_file(grn_file)
}

if(!interactive()){
  args <- getArgs()
  main(args)
}

