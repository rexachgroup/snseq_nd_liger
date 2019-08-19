# Damon Polioudakis
# 2019-02-09
# Run Seurat

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.6.0
#  module load python/3.7.2
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
# require(cellrangerRkit)
# require(loomR)
library(reticulate)
reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)
reticulate::py_config()   # check to make sure it is configured
source("function_library.R")
source("ggplot_theme.R")
# require(xlsx)

## Inputs
in_10x <- "/u/project/geschwind/drewse/g_singlecell/cellranger/data/20190624/aggr_20190625/outs/filtered_feature_bc_matrix/"
# metadata
p1_mt_df <- read_csv("../metadata/PGC_group1_071018.csv")
p2_mt_df <- read_csv("../metadata/PCG_round2.csv")
p3_mt_df <- read_csv("../metadata/PCG_round3_110918.csv")
p4_mt_df <- read_csv("../metadata/PCG_round4samples.csv")
# biomart gene info
bm_tb <- read_csv("../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv") %>% as_tibble %>% select(-X1)

## Variables
script_name <- "seurat.R"
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
out_graph <- paste0("../analysis/seurat/", date, "/graphs/p1_p2_p3_p4_seurat_")
out_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/", date
  , "/cellranger_v3/p1_p2_p3_p4_filtered_rmp27.rdat")
out_seurat_test <- gsub(".rdat", "_test.rdat", out_seurat)
out_seurat_raw <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/", date
  , "/cellranger_v3/p1_p2_p3_p4_raw.rdat")

# Make directories
# dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_seurat), recursive = TRUE)
################################################################################

### Main function

main_function <- function(){

  p1_p2_p3_p4_raw_so <- make_seurat_object(
    downsample_data = FALSE, ds_n_genes = 6000, ds_n_cells = 3000)
  save(p1_p2_p3_p4_raw_so, file = out_seurat_raw)

  p1_p2_p3_p4_so <- qc_filter_genes_and_cells(
      seurat_obj = p1_p2_p3_p4_raw_so
      , libraries_to_remove = "P2_7"
      , min_number_genes = 200
      , max_percent_mito = 5
      , min_cells_per_gene = 3)

  p1_p2_p3_p4_so <- run_seurat_pipeline(p1_p2_p3_p4_so)
  save(p1_p2_p3_p4_so, file = out_seurat)

  top_expressed_genes_tb <- find_top_cluster_expressed_genes(
    seurat_obj = p1_p2_p3_p4_so)
  save(p1_p2_p3_p4_so, top_expressed_genes_tb, file = out_seurat)
  p1_p2_p3_p4_so <- subset(p1_p2_p3_p4_so, downsample = 250)
  save(p1_p2_p3_p4_so, top_expressed_genes_tb, file = out_seurat_test)

  # load(paste0(
  #   "/u/flashscratch/d/dpolioud/seurat/20190427"
  #   , "/p1_p2_p3_p4_filtered_rmp27.rdat"))
  #
  # cluster_enriched_tb <- find_cluster_enriched_genes(seurat_obj = p1_p2_p3_p4_so)
  # save(p1_p2_p3_p4_so, top_expressed_genes_tb, cluster_enriched_tb, file = out_seurat)

}
################################################################################

### make to seurat object from expression matrix and metadata

make_seurat_object <- function(
  downsample_data = FALSE, ds_n_genes = 6000, ds_n_cells = 3000){

  print("make_seurat_object")

  # input 10X data and create seurat object
  # down-sample cells and genes for testing
  if (downsample_data == TRUE){
    print("down-sampling data")
    seurat_obj <-
      Read10X(data.dir = in_10x) %>%
        .[unique(c(
            sample(1:nrow(.), ds_n_genes)
            # include MT genes
            , grep(pattern = "^MT-", x = rownames(.))))
          , sample(1:ncol(.), ds_n_cells)] %>%
        CreateSeuratObject(counts = ., project = "p1-5_c1-5")
  # run full dataset
  } else {
    expr_M <- Read10X(data.dir = in_10x)
    print(str(expr_M))
    seurat_obj <- CreateSeuratObject(counts = expr_M, project = "p1-5_c1-5")
    rm(expr_M)
  }

  ## Add metadata

  format_metadata <- function(metadata_tb, vars_to_rename = NULL){
    metadata_tb <- metadata_tb %>%
      {if(! is.null(vars_to_rename)){rename(., !!vars_to_rename)} else .} %>%
      mutate("Sample #" = as.character(.$"Sample #")) %>%
      mutate(Prep = as.character(Prep)) %>%
      clean_variable_names()
    return(metadata_tb)
  }

  format_metadata_p1to5_c1to5 <- function(){
    print("format_metadata_p1to5_c1to5")
    # combine metadata from Jessica
    metadata_tb <- format_metadata(p1_mt_df) %>%
      bind_rows(., format_metadata(p2_mt_df)) %>%
      bind_rows(., format_metadata(p3_mt_df)) %>%
      bind_rows(., format_metadata(p4_mt_df)) %>%
      rename(., rin = "rin_after_staining_and_lcm")
    metadata_tb$region <- gsub("PreCG", "preCG", metadata_tb$region)

    tmp <- format_metadata(p5_c1to5_mt_df
        , vars_to_rename = c("Sample #" = "Sample.."
          , LIBRARY_ID = "Sample or Library Name"
          , Prep = "libraryBatch"
          # , Targeted_cellcount = "cell_counts"
          , targeted_cellcount = "targetcellcount(loaded)"
          , weight_g = "weight"
          ))
    tmp$targeted_cellcount <- gsub(",", "", tmp$targeted_cellcount) %>% as.numeric
    tmp$region <- gsub("PCG", "preCG", tmp$region)

    # clean columns
    metadata_tb <- bind_rows(metadata_tb, tmp) %>%
      select(-c(x14, x1, x19, x15, library_id_1, weightdissected, npdx1
        , cell_counts))

    # add cell ranger alignment metrics
    metadata_tb <- left_join(metadata_tb
      , cr_metrics_summary %>%
        # fix library ids
        # The core numbered P3 1 through 7 which matches in sequence to our 2 through 8
        mutate(ID = str_replace_all(ID, c(
          "P2_7B" = "P2_7_1B"
          , "C3_3_PC1_6" = "C3_3"
          , "P3_7" = "P3_8"
          , "P3_6" = "P3_7"
          , "P3_5" = "P3_6"
          , "P3_4" = "P3_5"
          , "P3_3" = "P3_4"
          , "P3_2" = "P3_3"
          , "P3_1" = "P3_2"))) %>%
        rename(library_id = "ID")
      , by = "library_id")

    # cell ranger adds numerical id tags to aggregated samples based on the order the samples were supplied in the aggregate csv
    # use the cell ranger numerical id tags and aggregate csv to match with metadata
    # first match numerical ids to aggregate csv and cleanup names in aggregate csv
    sample_id <- seurat_obj@assays$RNA@counts %>%
      colnames %>%
      gsub(".*-", "", .)
    cell_ranger_id_lib_id_key_tb <- left_join(
      sample_id %>% enframe()
      , cell_ranger_id_key_tb %>% mutate(value = rownames(cell_ranger_id_key_tb))
      ) %>%
        # fix library ids
        # The core numbered P3 1 through 7 which matches in sequence to our 2 through 8
        mutate(library_id = str_replace_all(library_id, c(
          "P2_7B" = "P2_7_1B"
          , "C3_3_PC1_6" = "C3_3"
          , "P3_7" = "P3_8"
          , "P3_6" = "P3_7"
          , "P3_5" = "P3_6"
          , "P3_4" = "P3_5"
          , "P3_3" = "P3_4"
          , "P3_2" = "P3_3"
          , "P3_1" = "P3_2"
      )))

    cell_ranger_id_lib_id_key_tb$library_id %>% unique
    # now match with metadata
    idx <- match(
      cell_ranger_id_lib_id_key_tb$library_id
      , metadata_tb$library_id)
    metadata_tb <- metadata_tb[idx, ]
    metadata_df <- metadata_tb %>% as.data.frame()
    metadata_df$cell_ids <- rownames(seurat_obj@meta.data)

    # Percent mito
    sum_mito <- grep(
        pattern = "^MT-", x = rownames(seurat_obj@assays$RNA@counts)) %>%
      seurat_obj@assays$RNA@counts[., ] %>%
      colSums()
    sum_all <- colSums(seurat_obj@assays$RNA@counts[ , ])
    percent_mito <- (sum_mito / sum_all) * 100
    metadata_df$percent_mito <- percent_mito

    # Number genes detected per cell
    metadata_df$number_genes <- seurat_obj@meta.data$nFeature_RNA

    # Number UMI per cell
    metadata_df$number_umi <- seurat_obj@meta.data$nCount_RNA

    # Cell counts per library id before filtering
    raw_cell_counts_per_lib <- metadata_df$library_id %>% table %>% as_tibble %>% rename()
    metadata_df <- left_join(metadata_df, raw_cell_counts_per_lib, by = c("library_id" = ".")) %>% rename(raw_cell_counts_per_lib = "n")

    rownames(metadata_df) <- metadata_df$cell_ids

    # # mean CDS length
    # mean_cds_length <- raw_loom$map(
    #   FUN = function(x){
    #     print(class(x))
    #     (x * raw_loom[["row_attrs/cds_length"]][]) %>%
    #     rowSums(., na.rm = TRUE) %>% as_tibble
    #   }
    #   , MARGIN = 2, chunk.size = 5000, dataset.use = "matrix",
    #   display.progress = TRUE
    #   ) %>%
    #   mutate(mean_cds_length = value /
    #     raw_loom[["col_attrs/number_umi"]][]) %>%
    #   pull(mean_cds_length)
    # raw_loom$add.col.attribute(
    #   list(mean_cds_length = mean_cds_length))
    #
    # # mean gc content
    # mean_gc_content <- raw_loom$map(
    #   FUN = function(x){
    #     print(class(x))
    #     (x * raw_loom[["row_attrs/percentage_gc_content"]][]) %>%
    #     rowSums(., na.rm = TRUE) %>% as_tibble
    #   }
    #   , MARGIN = 2, chunk.size = 5000, dataset.use = "matrix",
    #   display.progress = TRUE
    #   ) %>%
    #   mutate(mean_gc_content = value /
    #     raw_loom[["col_attrs/number_umi"]][]) %>%
    #   pull(mean_gc_content)
    # raw_loom$add.col.attribute(
    #   list(mean_gc_content = mean_gc_content))

    return(metadata_df)
  }

  format_gene_metadata <- function(){
    # Number cells per gene
    metadata_df$number_cells_per_gene <- rowSums(seurat_obj@assays$RNA@counts > 0)
  }

  metadata_df <- format_metadata_p1to5_c1to5()
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata_df)

  return(seurat_obj)
}
################################################################################

### QC and filter cells

qc_filter_genes_and_cells <- function(
  seurat_obj
  , libraries_to_remove = NULL
  , min_number_genes = 200
  , max_percent_mito = 5
  , min_cells_per_gene = 0){

  print("qc_filter_genes_and_cells")

  cells_to_keep <- seurat_obj@meta.data %>%
    # rownames_to_column(var = "cell_ids") %>%
    as_tibble() %>%
    filter(
      number_genes > min_number_genes
      , percent_mito < max_percent_mito) %>%
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

  return(seurat_obj)
}
################################################################################

run_seurat_pipeline <- function(seurat_obj){

  print("run_seurat_pipeline")

  print("NormalizeData")
  seurat_obj <- NormalizeData(object = seurat_obj, display.progress = TRUE)

  print("FindVariableGenes")
  seurat_obj <- FindVariableFeatures(
    seurat_obj, selection.method = "mean.var.plot", display.progress = TRUE)

  print("ScaleData")
  seurat_obj <- ScaleData(seurat_obj
      , do.scale = TRUE, do.center = TRUE, features = rownames(seurat_obj)
      # , vars.to.regress = "library_id"
      , display.progress = TRUE, num.cores = 8, do.par = TRUE)

  print("RunPCA")
  seurat_obj <- RunPCA(seurat_obj, online.pca = TRUE, npcs = 100
      , do.print = TRUE, pcs.print = 1:5, genes.print = 5
      , display.progress = TRUE)

  print("RunTSNE")
  seurat_obj <- RunTSNE(
    seurat_obj, dims.use = 1:100, do.fast = TRUE, nthreads = 8
    , tsne.method = "Rtsne", reduction = "pca", max_iter = 2000)

  print("RunUMAP")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:100)

  print("FindNeighbors")
  seurat_obj <- FindNeighbors(object = seurat_obj
    , reduction = "pca", dims = 1:100, nn.eps = 0, k.param = 30)

  print("FindClusters")
  for(i in c(0.4,0.5,0.6,0.7,0.8)){
    print(paste0("clustering for resolution: ", i))
    seurat_obj <- FindClusters(
      object = seurat_obj, resolution = i, n.start = 100)
  }

  return(seurat_obj)
}
################################################################################

### cluster enriched genes

convert_gene_symbols_to_ensembl_ids <- function(
  gene_symbols
  , in_biomart_gene_info = "../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"){
  # Converts vector of gene symbols to ensembl IDs
  # gene symbols must be character format
  # Need to load:
  biomart_df <- read.csv(in_biomart_gene_info, header = TRUE)
  print("convert_gene_symbols_to_ensembl_ids")
  idx <- match(gene_symbols, biomart_df$hgnc_symbol)
  ens <- biomart_df$ensembl_gene_id[idx]
  gene_symbols[! is.na(ens)] <- as.character(ens[! is.na(ens)])
  return(ens)
}

# Filter expression matrix by:
# percent of cells in cluster gene is expressed in
# fold change of gene in cluster versus all other cells
filter_expression_matrix_for_de <- function(
  expr_m
  , min_percent = NULL
  , fold_change = NULL
  , cluster_id = NULL
  , cell_cluster_key
  , cell_id = NULL) {

  if (! is.null(min_percent)) {
    # Expressed > 0 counts in > X% of cells in cluster
    if (! is.null(cluster_id)) {
      # Subset expression matrix to cluster
      cdf <- expr_m[ ,cell_cluster_key == cluster_id]
    }
    if (! is.null(cell_id)) {
      # Subset expression matrix to cluster
      cdf <- expr_m[ ,colnames(expr_m) %in% cell_id]
    }
    # Expressed > 0 counts in > X% of cells in cluster
    idxp <- (rowSums(cdf > 0) / ncol(cdf)) > (min_percent / 100)
    print(paste0("Genes expressed in > ", min_percent, "% of cells in cluster"))
    print(table(idxp))
  } else {
    idxp <- rep(TRUE, nrow(expr_m))
  }

  ### not currently generalized outside of Seurat
  if (! is.null(fold_change)) {
    # Fold change > Y of gene in cluster versus all other cells
    if (! is.null(cluster_id)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,cell_cluster_key == cluster_id]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! cell_cluster_key == cluster_id]
    }
    if (! is.null(cell_id)) {
      # Subset expression matrix to cluster
      cdf <- noCentExM[ ,colnames(expr_m) %in% cell_id]
      # Subset expression matrix to all other cells
      ndf <- noCentExM[ , ! colnames(expr_m) %in% cell_id]
    }
    # Fold change
    v1 <- rowMeans(cdf) - rowMeans(ndf)
    idxf <- v1 > fold_change
    print(paste0("Genes > ", fold_change, " fold change in cluster versus all other cells"))
    print(table(idxf))
  } else {
    idxf <- rep(TRUE, nrow(expr_m))
  }

  # Filter expr_m
  expr_m <- expr_m[idxp & idxf, ]
  return(expr_m)
}

## Function: DE Linear model
# terms_df:
# ExpCondition RIN.y     Seq.PC1      Seq.PC2
# 1            CP   8.4  0.04792498 -0.090448567
# 2            CP   8.1  0.53502697 -0.287629654
# 3            CP   8.0  0.18824922 -0.155651102
# 4            VZ   8.4  0.02529722 -0.100858264
# 5            VZ   8.7  0.45139297  0.856908177
# 6            VZ   9.1  0.27861748 -0.248868277
# mod: "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2"
run_de_with_linear_model <- function(expr_m, terms_df, mod) {
  print("run_de_with_linear_model")
  lmmod <- apply(expr_m, 1
    , function(y) {
      mod <- as.formula(mod)
      lm(mod, data = terms_df)})
  coefmat <- matrix(NA, nrow = nrow(expr_m)
    , ncol = length(coef(lmmod[[1]])))
  pvalmat <- matrix(NA, nrow = nrow(expr_m)
    , ncol = length(summary(lmmod[[1]])[[4]][ ,4]))
  colnames(coefmat) <- names(coef(lmmod[[1]]))
  rownames(coefmat) <- rownames(expr_m)
  colnames(pvalmat) <- names(summary(lmmod[[1]])[[4]][ ,4])
  rownames(pvalmat) <- rownames(expr_m)
  for (i in 1:nrow(expr_m)) {
    if (i%%100 == 0) {cat(".")}
    coefmat[i, ] <- coef(lmmod[[i]])
    pvalmat[i, ] <- summary(lmmod[[i]])[[4]][ ,4]
  }
  lm_coef_pval_l <- list(coefmat = coefmat, pvalmat = pvalmat)
  return(lm_coef_pval_l)
}

# Format output of linear model into data frame
format_lm_de <- function(
  lm_coef_pval_l, expr_m, cluster_id, cell_cluster_key, cluster_annot = NULL) {
  print("format_lm_de")
  # Combine log2 fold changes, p-values
  de_df <- data.frame(gene = row.names(lm_coef_pval_l$coefmat)
    , log2_fold_change = lm_coef_pval_l$coefmat[ ,2]
    , pvalue = lm_coef_pval_l$pvalmat[ ,2])
  # Add cluster ID
  de_df$Cluster <- cluster_id
  # Add cluster annotation
  if (! is.null(cluster_annot)){
    de_df$Cluster_Annot <- cluster_annot
  }
  # Percent of cells in cluster expressing gene > 0 counts
  cdf <- expr_m[
    row.names(expr_m) %in% de_df$gene, cell_cluster_key == cluster_id]
  de_df$percent_cluster <- (rowSums(cdf > 0) / ncol(cdf)) * 100
  de_df$percent_cluster <- round(de_df$percent_cluster, 1)
  # Percent of all cells expressing gene > 0 counts
  de_df$percent_all <- (rowSums(expr_m[
    row.names(expr_m) %in% de_df$gene, ] > 0)
    / ncol(expr_m[row.names(expr_m) %in% de_df$gene, ])) * 100
  de_df$percent_all <- round(de_df$percent_all, 1)
  # Order by log fold change
  de_df <- de_df[order(-de_df$log2_fold_change), ]
  return(de_df)
}

run_de_with_lm <- function(cluster_id, seurat_obj){

  print("run_de_with_lm")
  print(paste0("Calculating DE for cluster: ", cluster_id))

  # Filter cells
  expr_m <- filter_expression_matrix_for_de(
    expr_m = GetAssayData(seurat_obj, slot = "data")
    , min_percent = 10
    # , fold_change = 0
    , cluster_id = cluster_id
    , cell_cluster_key = Idents(seurat_obj)
    )

  # DE Linear model
  terms_df <- data.frame(seurat_obj[[c("RNA_snn_res.0.6"), drop = FALSE]])
  terms_df <- terms_df[row.names(terms_df) %in% colnames(expr_m), , drop = FALSE]

  # Add term TRUE/FALSE cell is in cluster
  terms_df$cluster <- rep(FALSE, nrow(seurat_obj[[]]))
  terms_df$cluster[Idents(seurat_obj) == cluster_id] <- TRUE
  mod <- "y ~ cluster"
  lm_coef_pval_l <- run_de_with_linear_model(
    expr_m = expr_m, terms_df = terms_df, mod = mod)

  # Format LM output into data frame
  cluster_enriched_de_df <- format_lm_de(
    lm_coef_pval_l = lm_coef_pval_l
    , expr_m = expr_m
    , cluster_id = cluster_id
    , cell_cluster_key = Idents(seurat_obj))

  # Add ensembl
  cluster_enriched_de_df$ensembl <- convert_gene_symbols_to_ensembl_ids(
    gene_symbols = cluster_enriched_de_df$gene)

  # fdr_pvalue correct
  cluster_enriched_de_df$fdr_pvalue <- p.adjust(
    cluster_enriched_de_df$pvalue, method = "BH")
  # Check
  table(cluster_enriched_de_df$pvalue < 0.05)
  table(cluster_enriched_de_df$fdr_pvalue < 0.05)

  # Format
  cluster_enriched_de_df$cluster <- cluster_id
  # Order columns
  cluster_enriched_de_df <- cluster_enriched_de_df[ ,c(
    "cluster", "gene", "ensembl", "log2_fold_change", "pvalue"
    , "fdr_pvalue", "percent_cluster", "percent_all")]
  # convert factors to characters
  cluster_enriched_de_df <- cluster_enriched_de_df %>%
    mutate_if(is.factor, as.character)

  return(cluster_enriched_de_df)
}

find_cluster_enriched_genes <- function(seurat_obj){

  print("find_cluster_enriched_genes")

  # Run DE
  clusters <- Idents(seurat_obj) %>% unique() %>% sort()
  cluster_enriched_de_tb <- map(clusters, .f = function(cluster_id){
    run_de_with_lm(cluster_id = cluster_id, seurat_obj = seurat_obj)
  }) %>% bind_rows() %>% as_tibble()

  # cluster_enriched_df <- FindAllMarkers(
  #   seurat_obj, test.use = "negbinom", latent.vars = "library_id"
  #   , only.pos = TRUE)

  return(cluster_enriched_de_tb)
}
################################################################################

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
################################################################################

### Run main function

main_function()
################################################################################

### CCA

run_cca <- function(){

  print("run_cca")

  ## Convert Loom to Seurat object

  filtered_loom <- connect(filename = out_filtered_loom, mode = "r+")
  # raw_loom <- connect(filename = out_raw_loom, mode = "r+")

  p1_p2_p3_p4_so <- Convert(from = filtered_loom, to = "seurat"
    , raw.data = "matrix"
    , gene.names = "row_attrs/gene_names"
    , cell.names = "col_attrs/cell_names"
    , norm.data = "layers/norm_data"
    , scale.data = "layers/scale_data"
    , gene.means = "row_attrs/gene_means"
    , gene.dispersion = "row_attrs/gene_dispersion"
    , gene.scaled = "row_attrs/gene_dispersion_scaled"
    , var.genes = "row_attrs/var_genes"
    )

  filtered_loom$close_all()


  ## CCA

  lib_ids <- unique(p1_p2_p3_p4_so@meta.data$library_id)

  cell_ids_1 <- p1_p2_p3_p4_so@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble %>%
    filter(library_id == lib_ids[1]) %>%
    pull(cell_id)

  ss_1_so <- SubsetData(p1_p2_p3_p4_so, cells.use = cell_ids_1)

  for(i in 2:length(lib_ids)){
  # for(i in 2:6){
    lib_id <- lib_ids[[i]]
    print(lib_id)

    cell_ids_2 <- p1_p2_p3_p4_so@meta.data %>% rownames_to_column("cell_id") %>% as_tibble %>% filter(library_id %in% lib_id) %>% pull(cell_id)

    ss_2_so <- SubsetData(p1_p2_p3_p4_so, cells.use = cell_ids_2)

    ss_1_so <- RunCCA(object = ss_1_so, object2 = ss_2_so, num.cc = 100)

  }


  smartseq2 <- CreateSeuratObject(raw.data = smartseq2.data)
  smartseq2 <- FilterCells(smartseq2, subset.names = "nGene", low.thresholds = 2500)
  smartseq2 <- NormalizeData(smartseq2)
  smartseq2 <- FindVariableGenes(smartseq2, do.plot = F, display.progress = F)
  smartseq2 <- ScaleData(smartseq2)
  smartseq2@meta.data$tech <- "smartseq2"

  # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  ob.list <- list(celseq, celseq2, fluidigmc1, smartseq2)
  genes.use <- c()
  for (i in 1:length(ob.list)) {
    genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
  }
  genes.use <- names(which(table(genes.use) > 1))
  for (i in 1:length(ob.list)) {
    genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
  }

  # Run multi-set CCA
  pancreas.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)

  # CC Selection
  MetageneBicorPlot(pancreas.integrated, grouping.var = "tech", dims.eval = 1:15)

  # Run rare non-overlapping filtering
  pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
                                         grouping.var = "tech", dims.use = 1:10)
  pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
                                             accept.low = 0.5)

  # Alignment
  pancreas.integrated <- AlignSubspace(pancreas.integrated,
                                       reduction.type = "cca",
                                       grouping.var = "tech",
                                       dims.align = 1:10)

  # t-SNE and Clustering
  pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
                                      dims.use = 1:10, save.SNN = T, resolution = 0.4)
  pancreas.integrated <- RunTSNE(pancreas.integrated,
                                 reduction.use = "cca.aligned",
                                 dims.use = 1:10)

  # Visualization
  TSNEPlot(pancreas.integrated, do.label = T)






  lib_ids <- unique(p1_p2_p3_p4_so@meta.data$library_id)

  cell_ids_1 <- p1_p2_p3_p4_so@meta.data %>%
    rownames_to_column("cell_id") %>%
    as_tibble %>%
    filter(library_id == lib_ids[1]) %>%
    pull(cell_id)

  ss_1_so <- SubsetData(p1_p2_p3_p4_so, cells.use = cell_ids_1)

  for(i in 2:length(lib_ids)){
  # for(i in 2:6){
    lib_id <- lib_ids[[i]]
    print(lib_id)

    cell_ids_2 <- p1_p2_p3_p4_so@meta.data %>% rownames_to_column("cell_id") %>% as_tibble %>% filter(library_id %in% lib_id) %>% pull(cell_id)

    ss_2_so <- SubsetData(p1_p2_p3_p4_so, cells.use = cell_ids_2)

    ss_1_so <- RunCCA(object = ss_1_so, object2 = ss_2_so, num.cc = 100)

  }

  save(ss_1_so, file = out_seurat)

  ss_1_so <- AlignSubspace(ss_1_so, reduction.type = "cca"
    , grouping.var = "library_id", dims.align = 1:75)

  ss_1_so <- RunTSNE(ss_1_so, reduction.use = "cca.aligned"
    , dims.use = 1:75, do.fast = TRUE)
  ss_1_so <- FindClusters(ss_1_so, reduction.type = "cca.aligned"
    , resolution = 0.8, dims.use = 1:75)

  save(ss_1_so, file = out_seurat)

  ## Plot
  bicor_data <- MetageneBicorPlot(ss_1_so, grouping.var = "library_id"
    , dims.eval = 1:100, display.progress = TRUE, return.mat = TRUE)
  ggplot(bicor_data, aes(x = cc, y = abs(bicor))) +
    geom_smooth(aes(col = Group, linetype = Group), se = FALSE) +
    scale_linetype_manual(values = rep(1:6,100)) +
    ylab(paste0("Shared correlation strength")) +
    xlab("CC")
  ggsave(paste(out_graph, "cca_metagene_bicor_plot.png"))

  as_tibble(ss_1_so@meta.data) %>%
    mutate(tsne_1 = ss_1_so@dr$tsne@cell.embeddings[,1]) %>%
    mutate(tsne_2 = ss_1_so@dr$tsne@cell.embeddings[,2]) %>%
  with(., plot_tsne_colored_by_variable(
    tsne_1 = tsne_1, tsne_2 = tsne_2, variable_value = library_id
    , facet_variable = "library_id", title = "CCA aligned by library"))
  ggsave(paste(out_graph, "tsne.png"))
    # , title = NULL, guide_size = 4, legend_title = NULL
    # , alpha = 0.5, size = 0.1
    # , expression_color_gradient = FALSE)
}
################################################################################

### UMAPs

# require(reticulate)
# sys <- import("sys")
# sys$version
# use_python("/u/local/apps/python/3.6.1/bin/python3")
# py_config()
#
# pbmc <- RunUMAP(p1_p2_p3_p4_so, reduction.use = "pca", dims.use = 1:10)
################################################################################
