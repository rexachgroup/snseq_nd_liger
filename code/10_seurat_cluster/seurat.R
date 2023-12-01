# Damon Polioudakis
# 2019-02-09
# Run Seurat

# Must load modules:
#  module load R/3.6.0
#  module load python/3.7.2

# submit as job:
# qsub -N srt -l h_data=256G,h_rt=72:00:00,highp qsub_r_script.sh -p seurat.R
###########################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(irlba)
require(tidyverse)
require(future)
library(reticulate)
reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)
reticulate::py_config()   # check to make sure it is configured
source("function_library.R")
source("seurat_function_library.R")

## Inputs
in_raw_so <- "/u/project/geschwind/drewse/g_singlecell/nucseq-nd/data/original_cells_w_doublets_nucseq_nd_pci_seurat_object.rds"
in_doublet_finder_so <- "/u/project/geschwind/drewse/g_singlecell/nucseq-nd/data/nucseq_nd_pci_dubs-rm_seurat_object.rds"
# biomart gene info
bm_tb <- read_csv("../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv") %>% as_tibble %>% select(-X1)
# cell type marker genes
marker_genes_refined_tb <- read_csv(
  "../resources/20191028_cell_markers_refined.csv")

## Variables
script_name <- "seurat.R"
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
# processed seurat object
out_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/", date, "/pci_filtered_seurat_object.rdat")
# processed seurat object down-sampled to make a test dataset
out_seurat_test <- gsub(".rdat", "_test.rdat", out_seurat)
# raw dataset as seurat object
out_seurat_raw <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/", date, "/pci_raw_seurat_object.rdat")
# raw dataset seurat object down-sampled to make a test dataset
out_seurat_raw_test <- gsub(".rdat", "_test.rdat", out_seurat_raw)

# Make directories
dir.create(dirname(out_seurat), recursive = TRUE)
###########################################################################

### Main function

main_function <- function(){

  print("main_function()")

  # load seurat object with doublets removed by andrew
  doublet_finder_so <- readRDS(in_doublet_finder_so)

  # doublet_finder_so <- subset(doublet_finder_so, downsample = 20000)

  # check dataset
  doublet_finder_so[[]] %>% count(clinical_dx, region, library_id) %>% as.data.frame() %>% print()

  # clean metadata and add my labels for consistency with subsequent scripts
  doublet_finder_so <- clean_and_format_metadata(doublet_finder_so)

  # check dataset
  doublet_finder_so[[]] %>% count(clinical_dx, region, library_id) %>% as.data.frame() %>% print()

  # QC filter
  nd_so <- qc_filter_genes_and_cells(
      seurat_obj = doublet_finder_so,
      # libraries_to_remove = "P2_7",
      min_number_genes = 200,
      max_number_genes = TRUE,
      max_percent_mito = 8,
      min_cells_per_gene = 3)
  # check
  print(dim(doublet_finder_so))
  print(dim(nd_so))
  rm(doublet_finder_so)

  # run seurat
  nd_so <- run_seurat_pipeline(nd_so)
  save(nd_so, file = out_seurat)

  # assign cell type to cells and clusters
  nd_so <- assign_cell_type_to_cell(
    seurat_obj = nd_so, marker_genes_tb = marker_genes_refined_tb)
  nd_so <- assign_cell_type_cluster(
    seurat_obj = nd_so, cluster_col_name = "RNA_snn_res.5")
  save(nd_so, file = out_seurat)

  ## downsample to make test datasets
  # downsample filtered dataset
  nd_so <- subset(nd_so, downsample = 250)
  # load raw dataset to downsample as well
  nd_raw_so <- readRDS(in_raw_so)
  # add my metadata labels for consistency with subsequent scripts
  nd_raw_so <- clean_and_format_metadata(nd_raw_so)
  # determine down-sample index for raw dataset
  cell_ids_raw_test_set <- c(
      colnames(nd_so),
      colnames(nd_raw_so)[sample(1:length(colnames(nd_raw_so)), 5000)]) %>%
    unique()
  # save full raw dataset
  save(nd_raw_so, file = out_seurat_raw)
  # down-sampled raw dataset
  nd_raw_so <- subset(nd_raw_so, cells = cell_ids_raw_test_set)
  # save down-sampled filtered dataset
  save(nd_so, file = out_seurat_test)
  # save down-sampled raw dataset
  save(nd_raw_so, file = out_seurat_raw_test)

  print("output paths:")
  print(out_seurat)
  print(out_seurat_test)
  print(out_seurat_raw)
  print(out_seurat_raw_test)

  print("main_function()")

}
###########################################################################

clean_and_format_metadata <- function(seurat_obj){

  print("clean_and_format_metadata()")

  # clean metadata and add my labels for consistency with subsequent scripts
  seurat_obj$number_genes <- seurat_obj[["nFeature_RNA"]]
  seurat_obj$number_umi <- seurat_obj[["nCount_RNA"]]
  seurat_obj$log_number_umi <- log(seurat_obj[["nCount_RNA"]])
  seurat_obj$percent_mito <- seurat_obj[["percent.mt"]]
  seurat_obj$cell_ids <- colnames(seurat_obj)
  seurat_obj$region <- seurat_obj[["region"]] %>%
    mutate(region = if_else(
      region == "PreCG", "preCG", if_else(
        region == "midInsula", "insula", region
    ))) %>% pull(region)

  print("end of... clean_and_format_metadata()")
  return(seurat_obj)
}
###########################################################################

### QC and filter cells

qc_filter_genes_and_cells <- function(
  seurat_obj,
  libraries_to_remove = NULL,
  min_number_genes = 200,
  max_number_genes = TRUE,
  max_percent_mito = 5,
  min_cells_per_gene = 0){

  print("qc_filter_genes_and_cells()")

  # set max_number_genes filter: mean + 3 SD
  if(max_number_genes == TRUE){
    mean_number_genes <- mean(seurat_obj[["number_genes"]]$number_genes)
    sd_number_genes <- sd(seurat_obj[["number_genes"]]$number_genes)
    max_number_genes <- mean_number_genes + (3 * sd_number_genes)
    print(paste0("max number genes filter (mean + 3 SD): ", max_number_genes))
  } else {
    max_number_genes <- max(seurat_obj[["number_genes"]]) + 1
  }

  # filter to cell ids to keep
  cells_to_keep <- seurat_obj@meta.data %>%
    # rownames_to_column(var = "cell_ids") %>%
    as_tibble() %>%
    filter(
      number_genes > min_number_genes,
      number_genes < max_number_genes,
      percent_mito < max_percent_mito) %>%
    # remove specific libraries
    {if(! is.null(libraries_to_remove)){
      filter(., ! library_id %in% libraries_to_remove)
    } else {.}} %>%
    select(cell_ids) %>%
    pull()

  genes_to_keep <- GetAssayData(seurat_obj, slot = "counts") %>%
    rowSums() %>%
    enframe(name = "gene", value = "number_cells") %>%
    filter(number_cells > min_cells_per_gene) %>%
    select(gene) %>%
    pull()

  seurat_obj <- subset(
    x = seurat_obj, cells = cells_to_keep, features = genes_to_keep)

  # Cell counts per library id after filtering
  cell_counts_per_lib <- seurat_obj[["library_id"]][,1] %>%
    table %>%
    as_tibble %>%
    rename(cell_counts_per_lib = "n")
  seurat_obj[["cell_counts_per_lib"]] <- left_join(
    seurat_obj[["library_id"]], cell_counts_per_lib
      , by = c("library_id" = ".")) %>%
    pull(cell_counts_per_lib)

  print("end of... qc_filter_genes_and_cells()")

  return(seurat_obj)
}
###########################################################################

run_seurat_pipeline <- function(seurat_obj){

  print("run_seurat_pipeline")

  print("NormalizeData")
  seurat_obj <- NormalizeData(object = seurat_obj, verbose = TRUE)

  print("FindVariableGenes")
  seurat_obj <- FindVariableFeatures(
    seurat_obj, selection.method = "vst", verbose = TRUE)
    # seurat_obj, selection.method = "mean.var.plot", verbose = TRUE)

  print("ScaleData")
  seurat_obj <- ScaleData(
    seurat_obj,
    do.scale = TRUE, do.center = TRUE,
    # features = rownames(seurat_obj),
    vars.to.regress = "log_number_umi",
    verbose = TRUE)

  print("RunPCA")
  seurat_obj <- RunPCA(seurat_obj, online.pca = TRUE, npcs = 200,
    do.print = TRUE, pcs.print = 1:5, genes.print = 5, verbose = TRUE)

  print("RunTSNE")
  seurat_obj <- RunTSNE(
    seurat_obj, dims.use = 1:100, do.fast = TRUE, nthreads = 8,
      tsne.method = "Rtsne", reduction = "pca", max_iter = 2000)

  print("RunUMAP")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:100)

  print("FindNeighbors")
  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "pca",
    dims = 1:100, nn.eps = 0, k.param = 30)

  print("FindClusters")
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = c(0.5, 1, 2, 3, 4, 5, 6),
    n.start = 100)

  return(seurat_obj)
}
###########################################################################

### top expressed genes by cluster

find_top_cluster_expressed_genes <- function(seurat_obj){

  print("find_top_cluster_expressed_genes")

  clusters <- sort(unique(Idents(seurat_obj)))

  top_expressed_genes_tb <-
    map(clusters, function(cluster){
      GetAssayData(seurat_obj, slot = "data")[
        , Idents(seurat_obj) == cluster] %>%
      rowMeans() %>%
      sort(decreasing = TRUE) %>%
      enframe() %>%
      as_tibble() %>%
      rename(gene = name, mean_expression = value) %>%
      mutate(cluster = cluster)
    }) %>%
    bind_rows()

  return(top_expressed_genes_tb)
}
###########################################################################

### Run main function

main_function()
###########################################################################
