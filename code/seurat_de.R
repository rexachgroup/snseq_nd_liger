# Damon Polioudakis
# 2019-05-29
# Determine enriched genes for a Seurat cluster

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

## Command args to input cluster ID
args <- commandArgs(trailingOnly = TRUE)
print(args)
cluster_id <- (as.numeric(args[1]) - 1)
# cluster_id <- 10
print(paste0("Cluster ID: ", cluster_id))

## Work paths
path_seurat <- args[2]
path_tmp <- args[5]

## Variables
script_name <- "seurat_de.R"
date <- format(Sys.Date(), "%Y%m%d")

seurat_obj <- as.character(args[3])
cluster_col_name <- as.character(args[4])
if (! is.na(args[7])) {
  compile_flag <- args[7]
} else {
  compile_flag <- FALSE
}

## Outputs

out_table <- args[6]

# Make directories
dir.create(dirname(path_tmp), recursive = TRUE)
################################################################################

### Main function

main_function <- function(){

  load(path_seurat)

  cluster_enriched_tb <- find_cluster_enriched_genes(
    seurat_obj = get(seurat_obj)
    , cluster_id = cluster_id
    , cluster_col_name = cluster_col_name)

  write.csv(cluster_enriched_tb, file = path_tmp, row.names = FALSE, quote = FALSE)

}

main_compile_function <- function(){

  load(path_seurat)

  cluster_enriched_tb <- map(
    list.files(dirname(path_tmp), full.names = TRUE), read.csv) %>% bind_rows()

  dir.create(dirname(out_table), recursive = TRUE)

  write.csv(cluster_enriched_tb, file = out_table, row.names = FALSE, quote = FALSE)

  save(top_expressed_genes_tb, cluster_enriched_tb, list = seurat_obj
    , file = path_seurat)
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

  # randomly subset cells not in cluster of interest to equal number in cluster of interest
  col_idxs <- c(0:ncol(expr_m))
  cell_keep_idx_cluster <- col_idxs[cell_cluster_key == cluster_id]
  # cell_keep_idx_noncluster <- sample(
  #   col_idxs[cell_cluster_key != cluster_id]
  #   , sum(cell_cluster_key == cluster_id))
  cell_keep_idx_noncluster <- sample(
    col_idxs[cell_cluster_key != cluster_id]
    , 10000)
  cell_keep_idx <- c(cell_keep_idx_cluster, cell_keep_idx_noncluster)

  # Filter expr_m
  expr_m <- expr_m[idxp & idxf, cell_keep_idx]

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

run_de_with_lm <- function(cluster_id, seurat_obj, cluster_col_name){

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
  terms_df <- data.frame(seurat_obj[[c(cluster_col_name), drop = FALSE]])
  # Add term TRUE/FALSE cell is in cluster
  terms_df$cluster <- rep(FALSE, nrow(seurat_obj[[]]))
  terms_df$cluster[Idents(seurat_obj) == cluster_id] <- TRUE
  terms_df <- terms_df[row.names(terms_df) %in% colnames(expr_m), , drop = FALSE]
  mod <- "y ~ cluster"
  lm_coef_pval_l <- run_de_with_linear_model(
    expr_m = expr_m, terms_df = terms_df, mod = mod)

  # Format LM output into data frame
  cell_cluster_key <- Idents(seurat_obj)[
    names(Idents(seurat_obj)) %in% colnames(expr_m)]
  cluster_enriched_de_df <- format_lm_de(
    lm_coef_pval_l = lm_coef_pval_l
    , expr_m = expr_m
    , cluster_id = cluster_id
    , cell_cluster_key = cell_cluster_key
  )

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

find_cluster_enriched_genes <- function(
  seurat_obj, cluster_id, cluster_col_name){

  print("find_cluster_enriched_genes")

  # Run DE
  cluster_enriched_de_tb <- run_de_with_lm(
      cluster_id = cluster_id
      , seurat_obj = seurat_obj
      , cluster_col_name = cluster_col_name
    ) %>% as_tibble()

  return(cluster_enriched_de_tb)
}
################################################################################

### Run main function

if (compile_flag == "TRUE") {
  main_compile_function()
} else {
  main_function()
}
################################################################################
