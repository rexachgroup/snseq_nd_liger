# Damon Polioudakis
# 2018-09-26
# Run Seurat on P1 and P2

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(irlba)
require(cowplot)
require(viridis)
require(tidyverse)
require(cellrangerRkit)
require(loomR)
source("function_library.R")
source("ggplot_theme.R")
# require(xlsx)

## Inputs
in_10x <- "../data/20181011/P1_P2/P1_P2/outs/filtered_gene_bc_matrices_mex/GRCh38-1.2.0_premrna/"
# metadata
p1_mt_df <- read_csv("../metadata/PGC_group1_071018.csv")
p2_mt_df <- read_csv("../metadata/PCG_round2.csv")
# biomart gene info
bm_tb <- read_csv("../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv") %>% as_tibble %>% select(-X1)

## Variables
script_name <- "p1_p2_seurat_loom.R"
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
out_graph <- paste0("../analysis/p1_p2_seurat_loom/graphs/", date, "/p1_p2_seurat_loom_")
out_raw_loom <- paste0(
  "/u/flashscratch/d/dpolioud/p1_p2_seurat_loom/", date
  , "/p1_p2_raw_p27removed_reglib.loom")
out_filtered_loom <- paste0(
  "/u/flashscratch/d/dpolioud/p1_p2_seurat_loom/", date
  , "/p1_p2_filtered_p27removed_reglib.loom")
out_filtered_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/p1_p2_seurat_loom/", date
  , "/p1_p2_filtered_p27removed_reglib.rdat")

# Make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_raw_loom), recursive = TRUE)
################################################################################

### Functions

main_function <- function(){
  make_loomr_object(
    downsample_data = FALSE, ds_n_genes = 6000, ds_n_cells = 3000)
  qc_filter_genes_and_cells(
    libraries_to_remove = "P2_7"
    , number_genes = 200
    , percent_mito = 5)
  read_depth_normalize()
  find_variable_genes()
  # center_scale_and_regress(vars_to_regress = c("library_id"))
  center_scale_and_regress()
  calculate_summary_stats()
  run_pca()
  run_tsne_and_clustering()
}
################################################################################

### make to loomR object from expression matrix and metadata

make_loomr_object <- function(
  downsample_data = FALSE, ds_n_genes = 6000, ds_n_cells = 3000){

  print("make_loomr_object")

  # input 10X data and create loomR
  # down-sample cells and genes for testing
  if (downsample_data == TRUE){
    print("down-sampling data")
    Read10X(data.dir = in_10x) %>%
      .[unique(c(
          sample(1:nrow(.), ds_n_genes)
          # include MT genes
          , grep(pattern = "^MT-", x = rownames(.))))
        , sample(1:ncol(.), ds_n_cells)] %>%
      create(out_raw_loom, .
        , overwrite = TRUE, display.progress = TRUE, do.transpose = TRUE
        , max.size = "4gb")
  # run full dataset
  } else {
    Read10X(data.dir = in_10x)[,] %>%
      create(out_raw_loom, .
        , overwrite = TRUE, display.progress = TRUE, do.transpose = TRUE
        , max.size = "4gb")
  }

  ## Add metadata

  raw_loom <- connect(filename = out_raw_loom, mode = "r+")

  clean_variable_names <- function(data){
    cleaned <- data %>%
      rename_all(
        funs(
          gsub("* ", "_", .) %>%
          gsub("\\.", "_", .) %>%
          gsub("\\(", "_", .) %>%
          gsub("\\)", "_", .) %>%
          gsub("\\+", "and", .) %>%
          gsub("#", "_number", .) %>%
          gsub("_$", "", .) %>%
          gsub("__", "_", .) %>%
          tolower
        )
      )
    return(cleaned)
  }

  format_metadata <- function(metadata_tb){
    metadata_tb <- metadata_tb %>%
        mutate("Sample #" = as.character(.$"Sample #")) %>%
        mutate(Prep = as.character(Prep)) %>%
        clean_variable_names()
    return(metadata_tb)
  }
  metadata <- format_metadata(p1_mt_df) %>%
    bind_rows(., format_metadata(p2_mt_df)) %>%
    # add column of cell ranger ids that are appended to cell_names
    # seem to be in the order as samples in the cell ranger aggregate args file
    mutate(cell_ranger_id = c(1:nrow(.)))
  sample_id <- gsub(".*-", "", raw_loom[["col_attrs/cell_names"]][])
  idx <- match(sample_id, metadata$cell_ranger_id)
  metadata <- metadata[idx, ]
  raw_loom$add.meta.data(metadata, overwrite = TRUE)

  # gene info
  row_dat_tb <- bm_tb %>%
    right_join(., raw_loom$get.attribute.df(MARGIN = 1) %>% as_tibble
      , by = c("hgnc_symbol" = "gene_names")) %>%
    distinct(., hgnc_symbol, .keep_all = TRUE)
  raw_loom$add.row.attribute(row_dat_tb, overwrite = TRUE)

  # Percent mito
  sum_mito <- grep(pattern = "^MT-", x = raw_loom$row.attrs$gene_names[]) %>%
    raw_loom[["matrix"]][,.] %>%
    rowSums()
  sum_all <- rowSums(raw_loom[["matrix"]][ , ])
  percent_mito <- (sum_mito / sum_all) * 100
  # Use add.row.attribute to add the IDs Note that if you want to overwrite an
  # existing value, set overwrite = TRUE
  raw_loom$add.col.attribute(list(percent_mito = percent_mito))

  # Number genes detected per cell
  number_genes <- rowSums(raw_loom[["matrix"]][ , ] > 0)
  raw_loom$add.col.attribute(list(number_genes = number_genes))

  # Number UMI per cell
  number_umi <- rowSums(raw_loom[["matrix"]][ , ])
  raw_loom$add.col.attribute(list(number_umi = number_umi))

  # Number cells per cell
  number_cells_per_gene <- colSums(raw_loom[["matrix"]][ , ])
  raw_loom$add.row.attribute(
    list(number_cells_per_gene = number_cells_per_gene))

  # mean CDS length
  mean_cds_length <- raw_loom$map(
    FUN = function(x){
      print(class(x))
      (x * raw_loom[["row_attrs/cds_length"]][]) %>%
      rowSums(., na.rm = TRUE) %>% as_tibble
    }
    , MARGIN = 2, chunk.size = 5000, dataset.use = "matrix",
    display.progress = TRUE
    ) %>%
    mutate(mean_cds_length = value /
      raw_loom[["col_attrs/number_umi"]][]) %>%
    pull(mean_cds_length)
  raw_loom$add.col.attribute(
    list(mean_cds_length = mean_cds_length))

  # mean gc content
  mean_gc_content <- raw_loom$map(
    FUN = function(x){
      print(class(x))
      (x * raw_loom[["row_attrs/percentage_gc_content"]][]) %>%
      rowSums(., na.rm = TRUE) %>% as_tibble
    }
    , MARGIN = 2, chunk.size = 5000, dataset.use = "matrix",
    display.progress = TRUE
    ) %>%
    mutate(mean_gc_content = value /
      raw_loom[["col_attrs/number_umi"]][]) %>%
    pull(mean_gc_content)
  raw_loom$add.col.attribute(
    list(mean_gc_content = mean_gc_content))

  raw_loom$close_all()
}
################################################################################

### Filter

qc_filter_genes_and_cells <- function(
  libraries_to_remove = NULL
  , number_genes_filter = 200
  , percent_mito_filter = 5
  ){

  print("qc_filter_genes_and_cells")

  raw_loom <- connect(filename = out_raw_loom, mode = "r+")

  cells_to_keep <- raw_loom$get.attribute.df(MARGIN = 2) %>%
    as_tibble %>%
    filter(number_genes > number_genes_filter
      , percent_mito < percent_mito_filter) %>%
    {if(! is.null(libraries_to_remove)){
      filter(., ! library_id %in% libraries_to_remove)
    } else {.}} %>%
    select(cell_names) %>%
    pull

  genes_to_keep <- raw_loom$get.attribute.df(MARGIN = 1) %>% as_tibble %>%
    filter(number_cells_per_gene > 3) %>% select(gene_names) %>% pull

  raw_loom[["matrix"]][
    raw_loom[["col_attrs/cell_names"]][] %in% cells_to_keep, ] %>%
    t() %>%
    as_tibble %>%
    mutate(gene_names = raw_loom[["row_attrs/gene_names"]][]) %>%
    filter(gene_names %in% genes_to_keep) %>%
    select(-gene_names) %>%
    create(out_filtered_loom, .
      , overwrite = TRUE, display.progress = TRUE, do.transpose = TRUE
      , max.size = "4gb")

  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")

  metadata <- raw_loom$get.attribute.df(MARGIN = 2) %>% as_tibble %>%
    filter(cell_names %in% cells_to_keep)
  filtered_loom$add.meta.data(metadata, overwrite = TRUE)

  gene_metadata <- raw_loom$get.attribute.df(MARGIN = 1) %>% as_tibble %>%
    filter(gene_names %in% genes_to_keep)
  filtered_loom$add.row.attribute(attribute = gene_metadata, overwrite = TRUE)

  raw_loom$close_all()
  filtered_loom$close_all()
}
################################################################################

### Normalize

read_depth_normalize <- function(){
  print("read_depth_normalize")
  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  # Saves as norm_data layer
  NormalizeData(object = filtered_loom, overwrite = FALSE, dtype = h5types$float
    , display.progress = TRUE
    # , chunk.size = 20000
    , chunk.dims = c(22361,22361), chunk.size = 22361
  )
  filtered_loom$close_all()
}
################################################################################

### Speed tests

# system.time(
# NormalizeData(object = raw_loom, overwrite = FALSE, dtype = h5types$float
#   , display.progress = TRUE
#   # , chunk.size = 20000
#   , chunk.dims = c(22361,22361), chunk.size = 22361
# ))
# # user  system elapsed
# # 193.799  27.620 223.456
#
# filtered_loom$close_all()
# raw_loom$close_all()
#
#
#
# system.time(
# p1_so <- Read10X(data.dir = in_10x) %>%
#   CreateSeuratObject(raw.data = .
#     , min.cells = 3, min.genes = 0,
#     , normalization.method = "LogNormalize", scale.factor = 10000
#     , project = "P1" , do.scale = FALSE, do.center = FALSE
#     , display.progress = TRUE)
# )
#
# system.time(
# p1_so <- NormalizeData(p1_so
#   , assay.type = "RNA", normalization.method = "LogNormalize"
#   , scale.factor = 10000, display.progress = TRUE
# ))
# # user  system elapsed
# # 12.725  17.567  48.239
################################################################################

### Detection of variable genes across the single cells

find_variable_genes <- function(){
  print("find_variable_genes")
  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  FindVariableGenes(filtered_loom, x.low.cutoff = 0.1
    , x.high.cutoff = 8, y.cutoff = 1, y.high.cutoff = Inf, num.bin = 20
    , chunk.size = 5000, normalized.data = "layers/norm_data"
    , overwrite = TRUE, display.progress = TRUE)
  filtered_loom$close_all()
}
################################################################################

### Center and scale data, regress out unwanted sources of variation

# Seurat constructs linear models to predict gene expression based on
# user-defined variables. The scaled z-scored residuals of these models are
# stored in the scale.data slot, and are used for dimensionality reduction and
# clustering. # note that this overwrites @scale.data. Therefore, if you intend
# to use ScaleData, you can set do.scale=F and do.center=F in the original
# object to save some time.

center_scale_and_regress <- function(vars_to_regress = NULL){
  # vars_to_regress format: c("percent_mito", "sample")

  print("center_scale_and_regress")

  ## Regress covariates and center scale

  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  # ScaleData(filtered_loom, genes.use = "row_attrs/var_genes", name =
  ScaleData(filtered_loom
    , genes.use = c(1:length(filtered_loom[["row_attrs/gene_names"]][]))
    , name = "scale_data", do.scale = TRUE, do.center = TRUE
    # , chunk.size = 22361, chunk.dims = c(22361, 22361)
    , vars.to.regress = vars_to_regress
    # , vars.to.regress = "percent_mito"
    , display.progress = TRUE, overwrite = TRUE)


  # ## No center or scale
  #
  # # Make dataframe of covariates to regress out
  # covDF <- data.frame(nUMI = so@meta.data$nUMI
  #   , librarylab = so@meta.data$librarylab
  #   , individual = so@meta.data$individual)
  # # Regress out confounding variables
  # RegressCovariates <- function (exM, covDF) {
  #   exRegCovM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
  #   rownames(exRegCovM) <- rownames(exM)
  #   colnames(exRegCovM) <- colnames(exM)
  #   # ncol(covDF)+1 when condition has 2 levels
  #   coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
  #   for (i in 1:nrow(exM)) {
  #     if (i%%1000 == 0) {print(paste("Done for", i, "genes..."))}
  #     mod <- lm(as.numeric(exM[i, ]) ~ ., data = covDF)
  #     # The full data - the undesired covariates
  #     exRegCovM[i,] <- coef(mod)[1] + mod$residuals
  #     # lmmod1 <- lm(as.numeric(exM[i, ]) ~ condition + age + sex + pmi, data = covDF)
  #   }
  #   return(exRegCovM)
  # }
  #
  # noCentExM <- RegressCovariates(so@data, covDF)

  filtered_loom$close_all()
}
################################################################################

### Calculate summary stats

calculate_summary_stats <- function(){

  print("calculate_summary_stats")

  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")

  # add gene summary stats for norm_data
  row_metadata <- tibble(
    mean_expression = filtered_loom$map(
      FUN = colMeans
      , MARGIN = 1, chunk.size = 2000, dataset.use = "matrix"
      , display.progress = TRUE
      )
    , median_expression = filtered_loom$map(
      FUN = function(x){
        x %>% as_tibble() %>% summarise_all(median) %>%
        unlist(., use.names = FALSE)
      }
      , MARGIN = 1, chunk.size = 200, dataset.use = "matrix"
      , display.progress = TRUE
      )
    , stdev_expression = filtered_loom$map(
      FUN = function(x){
        x %>% as_tibble() %>% summarise_all(sd) %>%
        unlist(., use.names = FALSE)
      }
      , MARGIN = 1, chunk.size = 200, dataset.use = "matrix"
      , display.progress = TRUE
      )
  )
  filtered_loom[["row_attrs/norm_data"]] <- row_metadata

  filtered_loom$close_all()
}
################################################################################

### Perform linear dimensional reduction

# Perform PCA on the scaled data. By default, the genes in object\@var.genes are
# used as input, but can be defined using pc.genes. We have typically found that
# running dimensionality reduction on genes with high-dispersion can improve
# performance. However, with UMI data - particularly after using ScaleData, we
# often see that PCA returns similar (albeit slower) results when run on much
# larger subsets of genes, including the whole transcriptome.

# Run PCA with the IRLBA package (iteratively computes the top dimensions,
# dramatic increase in speed since we only use a fraction of the PCs anyways) if
# you see the warning "did not converge–results might be invalid!; try
# increasing maxit or fastpath=FALSE", try increasing maxit

run_pca <- function(){
  print("run_pca")
  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  RunPCA(filtered_loom, online.pca = FALSE, pcs.compute = 100
    , do.print = TRUE, pcs.print = 1:5, genes.print = 5
    , display.progress = TRUE)
  filtered_loom$close_all()
}

run_ica <- function(){
  print("run_ica")
  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  RunICA(filtered_loom, ics.compute = 100, ics.print = 1:5, genes.print = 5)
  filtered_loom$close_all()
}
################################################################################

### Run Non-linear dimensional reduction (tSNE)

run_tsne_and_clustering <- function(){

  print("run_tsne_and_clustering")

  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")

  run_tsne_loom <- function(pcs, loom_col_name){
    print(paste0("run_tsne_loom PCs: ", pcs[1], "-", pcs[length(pcs)]
      , "; loom path: ", loom_col_name))
    RunTSNE(object = filtered_loom
      , reduction.use = "pca", dims.use = pcs, tsne.method = "Rtsne"
      , max_iter = 2000, nthreads = 4, overwrite = TRUE)
    col_to_add_l <- list()
    col_to_add_l[[loom_col_name]] = filtered_loom[["col_attrs/tsne_cell_embeddings"]][,]
    filtered_loom$add.col.attribute(col_to_add_l, overwrite = TRUE)
  }
  run_tsne_loom(pcs = 1:25, loom_col_name = "tsne_cell_embeddings_pc1to25")
  run_tsne_loom(pcs = 1:50, loom_col_name = "tsne_cell_embeddings_pc1to50")
  run_tsne_loom(pcs = 1:75, loom_col_name = "tsne_cell_embeddings_pc1to75")
  run_tsne_loom(pcs = 1:100, loom_col_name = "tsne_cell_embeddings_pc1to100")

  # FindClusters(object = filtered_loom
  #   , reduction.type = "pca", dims.use = 1:25, resolution = 0.8
  #   , save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)

  run_clustering_loom <- function(pcs, loom_col_name){
    print(paste0("run_clustering_loom PCs: ", pcs[1],"-", tail(pcs,1)
      , " loom path: ", loom_col_name))
  # By default, we perform 100 random starts for clustering and select the result
  # with highest modularity. You can lower this through the n.start parameter to
  # reduce clustering time.
      FindClusters(object = filtered_loom
        , reduction.type = "pca", dims.use = pcs, resolution = 0.8
        , save.SNN = FALSE, n.start = 10, nn.eps = 0.5, print.output = FALSE
        , force.recalc = TRUE)
      col_to_add_l <- list()
      col_to_add_l[[loom_col_name]] <- filtered_loom[["col_attrs/res.0.8"]][] %>%
        as.character
      filtered_loom$add.col.attribute(col_to_add_l, overwrite = TRUE)
  }
  run_clustering_loom(pcs = 1:25, loom_col_name = "cluster_pc1to25_res08")
  run_clustering_loom(pcs = 1:50, loom_col_name = "cluster_pc1to50_res08")
  run_clustering_loom(pcs = 1:75, loom_col_name = "cluster_pc1to75_res08")
  run_clustering_loom(pcs = 1:100, loom_col_name = "cluster_pc1to100_res08")

  filtered_loom$close_all()
}

## Resolution test

################################################################################

### Run main function

main_function()
################################################################################

### UMAPs

# require(reticulate)
# sys <- import("sys")
# sys$version
# use_python("/u/local/apps/python/3.6.1/bin/python3")
# py_config()
#
# pbmc <- RunUMAP(p1_p2_so, reduction.use = "pca", dims.use = 1:10)
################################################################################
