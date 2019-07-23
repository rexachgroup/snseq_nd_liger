# Damon Polioudakis
# 2019-02-09
# Run Seurat on P1 and P2

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
#  module load python/3.7.2
################################################################################

rm(list = ls())
set.seed(27)

.libPaths("/u/home/d/dpolioud/R/x86_64-pc-linux-gnu-library/3.4/")
require(methods)
require(Seurat)
require(Matrix)
require(irlba)
require(cowplot)
require(viridis)
require(cellrangerRkit)
require(tidyverse)
library(reticulate)
reticulate::use_python("/u/local/apps/python/3.7.2/bin/python3", required=TRUE)
reticulate::py_config()   # check to make sure it is configured
source("function_library.R")
source("ggplot_theme.R")
# require(xlsx)
sessionInfo()

## Inputs
in_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/20190427"
  , "/p1_p2_p3_p4_filtered_rmp27.rdat")

## Variables
script_name <- "seurat_subcluster.R"
date <- format(Sys.Date(), "%Y%m%d")
# date <- "20190221"

## Outputs
out_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/seurat_subcluster/", date
  , "/seuratv3_integration/seurat_subcluster_astro.rdat")
out_seurat_test <- gsub(".rdat", "_test.rdat", out_seurat)
out_table <- paste0(
  "../analysis/seurat_subcluster/", date
  , "/tables/seuratv3_integration/seurat_subcluster_astro_")
print(out_seurat)
print(out_table)

# Make directories
dir.create(dirname(out_seurat), recursive = TRUE)
# dir.create(dirname(out_table), recursive = TRUE)
################################################################################

### Main function

main_function <- function(){

  load(in_seurat)
  Idents(p1_p2_p3_p4_so) <- p1_p2_p3_p4_so[["RNA_snn_res.0.6"]]

  subset_so <- subset(p1_p2_p3_p4_so, idents = c(3,14,22))
    # , subset = clinical_dx == "Control")
  rm(p1_p2_p3_p4_so)

  # subset_so <- readRDS(in_seurat)
  subset_so <- run_seurat_cca_pipeline(subset_so, dims_use = 1:10)
  save(subset_so, file = out_seurat)

  # stash clustering of full dataset for later
  # subset_so$cluster_pc1to100_res06 <- subset_so[["RNA_snn_res.0.6"]]
  # subset_so <- run_seurat_integration_pipeline(subset_so, dims_use = 1:30)
  # Idents(subset_so) <- subset_so[["integrated_snn_res.0.4"]]
  # save(subset_so, file = out_seurat)
  # subset_so <- subset(subset_so, downsample = 500)
  # save(subset_so, top_expressed_genes_tb, file = out_seurat_test)

  # stash clustering of full dataset for later
  # subset_so$cluster_pc1to100_res06 <- subset_so[["RNA_snn_res.0.6"]]
  # subset_so <- run_seurat_pipeline(subset_so)
  # save(subset_so, file = out_seurat)

  # set cluster identities to desired clustering
  Idents(subset_so) <- subset_so[["RNA_snn_res.0.4"]]
  top_expressed_genes_tb <- find_top_cluster_expressed_genes(subset_so)
  save(subset_so, top_expressed_genes_tb, file = out_seurat)
  subset_so <- subset(subset_so, downsample = 500)
  save(subset_so, top_expressed_genes_tb, file = out_seurat_test)

  #
  # cluster_enriched_de_tb <- find_cluster_enriched_genes(subset_so)
  # save(subset_so, top_expressed_genes_tb, cluster_enriched_de_tb
  #   , file = out_seurat)
  # write.csv(x = cluster_enriched_de_tb
  #   , file = paste0(out_table, "cluster_enriched_vs_all.csv")
  #   , quote = FALSE, row.names = FALSE)

}
################################################################################

run_seurat_pipeline <- function(seurat_obj){

  print("run_seurat_pipeline")

  # print("NormalizeData")
  # seurat_obj <- NormalizeData(object = seurat_obj, display.progress = TRUE)

  print("FindVariableGenes")
  seurat_obj <- FindVariableFeatures(
    # seurat_obj, selection.method = "mean.var.plot", display.progress = TRUE)
    seurat_obj, selection.method = "vst", display.progress = TRUE)
  print(
    paste0("Number of variable genes: ", length(VariableFeatures(seurat_obj))))

  print("ScaleData")
  seurat_obj <- ScaleData(seurat_obj
      , do.scale = TRUE, do.center = TRUE, features = rownames(seurat_obj)
      # , vars.to.regress = "rin_after_staining_and_lcm"
      # , vars.to.regress = "library_id"
      , display.progress = TRUE, num.cores = 8, do.par = TRUE)

  print("RunPCA")
  seurat_obj <- RunPCA(seurat_obj, online.pca = FALSE, npcs = 100
      , do.print = TRUE, pcs.print = 1:5, genes.print = 5
      , display.progress = TRUE)

  print("RunTSNE")
  seurat_obj <- RunTSNE(
    seurat_obj, dims.use = 1:25, do.fast = TRUE, nthreads = 8
    , tsne.method = "Rtsne", reduction = "pca", max_iter = 2000)
  print("RunUMAP")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:25)

  print("FindNeighbors")
  seurat_obj <- FindNeighbors(object = seurat_obj
    , reduction = "pca", dims = 1:25, nn.eps = 0, k.param = 30)

  print("FindClusters")
  for(i in c(0.4,0.5,0.6,0.7,0.8)){
    print(paste0("clustering for resolution: ", i))
    seurat_obj <- FindClusters(
      object = seurat_obj, resolution = i, n.start = 100)
  }

  return(seurat_obj)
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
      Matrix::rowMeans() %>%
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
run_de_with_lm_with_linear_model <- function(expr_m, terms_df, mod) {
  print("run_de_with_lm_with_linear_model")
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
  lm_coef_pval_l <- run_de_with_lm_with_linear_model(
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

### CCA

run_seurat_cca_pipeline <- function(seurat_obj, dims_use = 1:10){

  print("run_cca")

  expr_m <- GetAssayData(seurat_obj, "counts")
  metadata_df <- seurat_obj[[]]

  detach("package:Seurat", unload = TRUE)
  .libPaths()
  .libPaths("/u/home/d/dpolioud/R/x86_64-pc-linux-gnu-library/3.4/Seurat_2.3.4")
  .libPaths()
  require(Seurat)

  seurat_obj <- CreateSeuratObject(raw.data = expr_m, project = "p1_p2_p3_p4_astro")
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata_df)

  lib_ids <- seurat_obj@meta.data["library_id"] %>% unique() %>% pull()

  seurat_obj_l <- map(.x = lib_ids, .f = function(lib_id){

    print(lib_id)

    cell_ids <- seurat_obj@meta.data["library_id"] %>%
      as_tibble(rownames = "cell_ids") %>%
      filter(library_id == lib_id) %>%
      pull(cell_ids)

    seurat_obj <- SubsetData(seurat_obj, cells = cell_ids)

    print("NormalizeData")
    seurat_obj <- NormalizeData(object = seurat_obj, display.progress = TRUE)

    print("FindVariableGenes")
    seurat_obj <- FindVariableGenes(
      seurat_obj, selection.method = "mean.var.plot", display.progress = TRUE)

    print("ScaleData")
    seurat_obj <- ScaleData(seurat_obj
        , do.scale = TRUE, do.center = TRUE
        , display.progress = TRUE, num.cores = 8, do.par = TRUE)

    return(seurat_obj)

  })

  # Determine genes to use for CCA, must be highly variable in at least 5 datasets
  ob.list <- seurat_obj_l
  genes.use <- c()
  for (i in 1:length(ob.list)) {
    genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 2000))
  }
  genes.use <- names(which(table(genes.use) > 4))
  for (i in 1:length(ob.list)) {
    genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
  }
  # 3369
  length(genes.use)

  # Run multi-set CCA
  seurat_obj <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 25)

  # CC Selection
  # MetageneBicorPlot(seurat_obj, grouping.var = "library_id", dims.eval = 1:100)
  #   ggsave(paste0("metagenebicorplot.png"))

  # Run rare non-overlapping filtering
  seurat_obj <- CalcVarExpRatio(object = seurat_obj
    , reduction.type = "pca"
    , grouping.var = "library_id"
    , dims.use = dims_use)

  # seurat_obj <- SubsetData(
  #   seurat_obj
  #   , subset.name = "var.ratio.pca"
  #   , accept.low = 0.5)

  # Alignment
  print("AlignSubspace")
  seurat_obj <- AlignSubspace(
    seurat_obj
    , reduction.type = "cca"
    , grouping.var = "library_id"
    , dims.align = dims_use)

  print("FindClusters")
  # t-SNE and Clustering
  seurat_obj <- FindClusters(
    seurat_obj
    , reduction.type = "cca.aligned"
    , dims.use = dims_use
    , save.SNN = T
    , resolution = 0.6)

  print("RunTSNE")
  seurat_obj <- RunTSNE(
    seurat_obj
    , reduction.use = "cca.aligned"
    , dims.use = dims_use)

  print("RunUMAP")
  seurat_obj <- RunUMAP(
    seurat_obj
    , dims = dims_use
    , reduction = "cca.aligned")

  return(seurat_obj)
  # Visualization
  # TSNEPlot(pancreas.integrated, do.label = T)
  # ggsave(paste0("cca_tsne.png"))

  # # code for seurat v3.0
  # ob.list <- seurat_obj_l
  # genes.use <- c()
  # for (i in 1:length(ob.list)) {
  #   genes.use <- c(genes.use, VariableFeatures(ob.list[[i]]))
  # }
  # genes.use <- names(which(table(genes.use) > 1))
  # for (i in 1:length(ob.list)) {
  #   genes.use <- genes.use[genes.use %in% rownames(GetAssayData(object = ob.list[[i]], slot = "scale.data"))]
  # }

  # # code for seurat v3.0 CCA
  # seurat_obj <- ob.list[[1]]
  # for (i in 2:length(ob.list)){
  #   seurat_obj <- RunCCA(
  #     seurat_obj, ob.list[[i]], genes.use = genes.use, num.ccs = 15)
  #   seurat_obj <- FindVariableFeatures(
  #     seurat_obj, selection.method = "mean.var.plot", display.progress = TRUE)
  # }



  #
  #
  #
  # smartseq2 <- CreateSeuratObject(raw.data = smartseq2.data)
  # smartseq2 <- FilterCells(smartseq2, subset.names = "nGene", low.thresholds = 2500)
  # smartseq2 <- NormalizeData(smartseq2)
  # smartseq2 <- FindVariableGenes(smartseq2, do.plot = F, display.progress = F)
  # smartseq2 <- ScaleData(smartseq2)
  # smartseq2@meta.data$tech <- "smartseq2"
  #
  # # Determine genes to use for CCA, must be highly variable in at least 2 datasets
  # ob.list <- list(celseq, celseq2, fluidigmc1, smartseq2)
  # genes.use <- c()
  # for (i in 1:length(ob.list)) {
  #   genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
  # }
  # genes.use <- names(which(table(genes.use) > 1))
  # for (i in 1:length(ob.list)) {
  #   genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
  # }
  #
  # # Run multi-set CCA
  # pancreas.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)
  #
  # # CC Selection
  # MetageneBicorPlot(pancreas.integrated, grouping.var = "tech", dims.eval = 1:15)
  #
  # # Run rare non-overlapping filtering
  # pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
  #                                        grouping.var = "tech", dims.use = 1:10)
  # pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
  #                                            accept.low = 0.5)
  #
  # # Alignment
  # pancreas.integrated <- AlignSubspace(pancreas.integrated,
  #                                      reduction.type = "cca",
  #                                      grouping.var = "tech",
  #                                      dims.align = 1:10)
  #
  # # t-SNE and Clustering
  # pancreas.integrated <- FindClusters(pancreas.integrated, reduction.type = "cca.aligned",
  #                                     dims.use = 1:10, save.SNN = T, resolution = 0.4)
  # pancreas.integrated <- RunTSNE(pancreas.integrated,
  #                                reduction.use = "cca.aligned",
  #                                dims.use = 1:10)
  #
  # # Visualization
  # TSNEPlot(pancreas.integrated, do.label = T)

}
################################################################################

### Run Seurat integration pipeline

run_seurat_integration_pipeline <- function(seurat_obj, dims_use = 1:10){

  print("run_seurat_integration_pipeline")

  seurat_obj_l <- SplitObject(seurat_obj, split.by = "prep")
  seurat_obj_l <- SplitObject(seurat_obj, split.by = "library_id")

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  print("Normalize each batch individually")
  for (i in 1:length(seurat_obj_l)) {
      seurat_obj_l[[i]] <- NormalizeData(seurat_obj_l[[i]], verbose = FALSE)
      seurat_obj_l[[i]] <- FindVariableFeatures(
        seurat_obj_l[[i]], selection.method = "vst"
        , nfeatures = 2000, verbose = FALSE)
  }

  #lapply(seurat_obj_l, function(seurat_obj){GetAssayData(seurat_obj, slot = "data") %>% dim})

  print("Find integration anchors")
  anchors <- FindIntegrationAnchors(
    object.list = seurat_obj_l, dims = 1:30, k.filter = 50)

  integrated_seurat_obj <- IntegrateData(anchorset = anchors, dims = 1:30)

  DefaultAssay(integrated_seurat_obj) <- "integrated"

  print("Scale data")
  integrated_seurat_obj <- ScaleData(
    integrated_seurat_obj, verbose = FALSE
    , features = rownames(integrated_seurat_obj))

  print("Run PCA")
  integrated_seurat_obj <- RunPCA(
    integrated_seurat_obj, npcs = 30, verbose = FALSE)

  print("RunTSNE")
  integrated_seurat_obj <- RunTSNE(
    integrated_seurat_obj, dims.use = 1:30, do.fast = TRUE, nthreads = 8
    , tsne.method = "Rtsne", reduction = "pca", max_iter = 2000)

  print("RunUMAP")
  integrated_seurat_obj <- RunUMAP(integrated_seurat_obj, dims = 1:30)

  print("FindNeighbors")
  integrated_seurat_obj <- FindNeighbors(object = integrated_seurat_obj
    , reduction = "pca", dims = 1:30, nn.eps = 0, k.param = 30)

  print("FindClusters")
  for(i in c(0.4,0.5,0.6,0.7,0.8)){
    print(paste0("clustering for resolution: ", i))
    integrated_seurat_obj <- FindClusters(
      object = integrated_seurat_obj, resolution = i, n.start = 100)
  }

  return(integrated_seurat_obj)
}
################################################################################

### Run main function

main_function()
################################################################################
