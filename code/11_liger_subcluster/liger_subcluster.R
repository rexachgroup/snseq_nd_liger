# Damon Polioudakis
# 2020-01-20
# Sub-clustering with Liger (2nd round of clustering)

# Must load modules:
#  module load R/3.6.0
#  module load python/3.7.2

# Sample qsub:
# qsub -N liger_sc -t 1-30 -m n -tc 5 -l h_data=128G,h_rt=24:00:00,highp qsub_r_script.sh -p liger_subcluster.R
###########################################################################

rm(list = ls())
set.seed(27)

# require(devtools)
# devtools::install_github('MacoskoLab/liger')
# devtools::install_github('satijalab/seurat-wrappers')

require(Seurat)
require(SeuratWrappers)
require(liger)
require(tidyverse)
require(reticulate)
reticulate::use_condaenv("../../conda/", required=TRUE)
reticulate::import("umap")
reticulate::py_config()   # check to make sure it is configured
source("function_library.R")
source("ggplot_theme.R")
# require(xlsx)
sessionInfo()

## Command line arguments
command_args <- commandArgs(trailingOnly = TRUE)
print(command_args)

## Inputs
in_seurat <- "../..//analysis/seurat/20200816/pci_filtered_seurat_object.rdat"

## Variables
script_name <- "liger_subcluster.R"
graph_subtitle <- "P1-5 C1-5 I1-5, 8% MT filter, regress out log(number UMI)"

## Outputs
out_seurat_path <- paste0(
  "../analysis/seurat_lchen/liger_subcluster/liger_subcluster_")
# Make directories
dir.create(dirname(out_seurat_path), recursive = TRUE)
###########################################################################

### Main function

main_function <- function(){

  load(in_seurat)

  ## Define subset arguments
  #i <- as.numeric(command_args[1])
  cmd_region <- command_args[[1]]
  cmd_cell_type <- command_args[[2]]
  cell_subset_tb <- expand.grid(
    # clinical_dx = c("AD", "bvFTD", "PSP-S"),
    region = nd_so[["region"]] %>% pull %>% unique,
    cluster_cell_type = nd_so[["cluster_cell_type"]] %>% pull %>% unique)
  # need to make variables global for seurat subset() function
  # assign("dx_subset", c(cell_subset_tb$clinical_dx[i] %>% as.character(), "Control"), envir = .GlobalEnv)
  assign("dx_subset", c("AD", "bvFTD", "PSP-S", "Control"), envir = .GlobalEnv)
  assign("region_subset", cmd_region, envir = .GlobalEnv)
  assign("cell_type_subset", cmd_cell_type %>% as.character(), envir = .GlobalEnv)
  print(dx_subset)
  print(region_subset)
  print(cell_type_subset)

  ## Outputs
  out_seurat <- paste0(
    out_seurat_path,
    paste0(dx_subset, collapse = "_"), "_", region_subset, "_",
    cell_type_subset,".rdat")
  out_seurat_test <- gsub(".rdat", "_test.rdat", out_seurat)

  print(out_seurat)

  ## Subset
  # make sure Seurat:subset == [VARIABLE] is not a metadata column
  print("sub-setting")
  subset_meta <- nd_so@meta.data %>%
      as_tibble %>%
      filter(cluster_cell_type == cmd_cell_type, region == cmd_region)
  #   subset_so <- subset(nd_so, subset = cluster_cell_type == cell_type_subset)
  #   subset_so <- subset(subset_so, subset = region == region_subset)
  #   subset_so <- subset(subset_so, subset = clinical_dx %in% dx_subset)
  subset_so <- subset(nd_so, cells = subset_meta$cell_ids)
  print("Full dataset:")
  FetchData(nd_so, slot = "data",
    vars = c("cluster_cell_type", "clinical_dx", "region")) %>%
    table() %>%
    print()
  print("Subset dataset:")
  FetchData(subset_so, slot = "data",
    vars = c("cluster_cell_type", "clinical_dx", "region")) %>%
    table() %>%
    print()
  rm(nd_so)

  ## Run liger suggestK and suggestLamba functions
  # run_liger_suggest_k_and_lambda(subset_so)

  ## Run liger pipeline
  # stash clustering of full dataset for later
  subset_so$cluster_full_datatset <- Idents(subset_so)
  liger_so <- run_liger_pipeline(subset_so, k = 20, lambda = 0.5)
  save(liger_so, file = out_seurat)

  ## Find top expressed genes
  # set cluster identities to desired clustering
  Idents(liger_so) <- liger_so[["liger_clusters"]]
  top_expressed_genes_tb <- find_top_cluster_expressed_genes(liger_so)
  save(liger_so, top_expressed_genes_tb, file = out_seurat)

  ## Output csv of sub-clusters, cell ids, and liger rdat file name
  # sub-clusters and cell ids
  subclusters_tb <- liger_so[[c("liger_clusters")]] %>%
    as_tibble(rownames = "cell_ids") %>%
    mutate(liger_rdat_file = basename(out_seurat))
  out_dir <- paste0(dirname(out_seurat), "/tables")
  out_file_name <- basename(out_seurat) %>%
    gsub(".rdat", paste0("_subcluster_ids.csv"), .)
  out_subcluster_ids_csv <- paste0(out_dir, "/", out_file_name)
  dir.create(out_dir, recursive = TRUE)
  print("Saving csv too:")
  print(out_subcluster_ids_csv)
  write_csv(subclusters_tb, path = out_subcluster_ids_csv)

  ## Make test dataset by down-sampling
  liger_so <- subset(liger_so, downsample = 500)
  save(liger_so, top_expressed_genes_tb, file = out_seurat_test)

  #
  # ## Determine cluster enriched genes
  # cluster_enrich_seurat_df <- FindAllMarkers(liger_so, test.use = "negbinom",
  #   latent.vars = "number_umi")
  # cluster_enrich_seurat_tb <- cluster_enrich_seurat_df %>% as_tibble()
  # save(liger_so, top_expressed_genes_tb, cluster_enrich_seurat_tb,
  #   file = out_seurat)

  print("end of main_function")

}
###########################################################################

### Run Liger suggest k and lambda

run_liger_suggest_k_and_lambda <- function(subset_so){

  print("run_liger_suggest_k_and_lambda")

  out_graph <- assign(
    "out_graph",
    paste0(
      dirname(out_seurat_path),
      paste0(dx_subset, collapse = "_"), "_", region_subset, "_",
      cell_type_subset, "_"),
    envir = .GlobalEnv)

  dir.create(dirname(out_graph), recursive = TRUE)

  subset_so <- NormalizeData(subset_so)
  subset_so <- FindVariableFeatures(subset_so)
  subset_so <- ScaleData(subset_so, split.by = "prep", do.center = FALSE,
    vars.to.regress = "number_umi")

  ligerex <- seuratToLiger(subset_so,
    combined.seurat = T, names = "use-projects",
    meta.var = "prep", assays.use = NULL, raw.assay = "RNA",
    remove.missing = T, renormalize = T, use.seurat.genes = T,
    num.hvg.info = NULL, use.idents = T, use.tsne = T, cca.to.H = F)

  ligerex <- scaleNotCenter(ligerex)
  suggestK(ligerex, k.test = seq(5, 100, 5)) # plot entropy metric to find an elbow that can be used to select the number of factors
  ggsave(paste0(out_graph, "suggest_k.png"))

  # The default value of lambda=5.0 usually provides reasonable results for most analyses, although the suggestLambda function can be used to determine a more appropriate lambda value for the desired level of dataset alignment.
  suggestLambda(ligerex, k = 50) # plot alignment metric to find an elbow that can be used to select the value of lambda
  ggsave(paste0(out_graph, "suggest_lambda.png"))
}
###########################################################################

### Run Liger pipeline

run_liger_pipeline <- function(subset_so, k = 50, lambda = 3){

  print("run_liger_pipeline")

  print(paste0("k = ", k))
  print(paste0("lambda = ", lambda))

  ligerex <- seuratToLiger(subset_so,
    combined.seurat = T, names = "use-meta",
    meta.var = "prep", assays.use = NULL, raw.assay = "RNA",
    remove.missing = F, renormalize = F, use.seurat.genes = F,
    num.hvg.info = NULL, use.idents = F, use.tsne = T, cca.to.H = F)

  ligerex <- normalize(ligerex)
  # select genes
  ligerex <- selectGenes(ligerex)
  ligerex <- scaleNotCenter(ligerex)

  # Take the lowest objective of 5 factorizations with different initializations
  # Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
  # The default value of lambda=5.0 usually provides reasonable results for most analyses, although the suggestLambda function can be used to determine a more appropriate lambda value for the desired level of dataset alignment.
  ligerex <- optimizeALS(ligerex, k = k, lambda = lambda, thresh = 5e-5, nrep = 3)
  # run quantile normalization and alignment
  ligerex <- quantileAlignSNF(ligerex) #SNF clustering and quantile alignment
  ligerex <- runUMAP(ligerex)

  # Save cells ordering from seurat
  seurat.cell.ordering <- rownames(subset_so@reductions$tsne@cell.embeddings)

  # Add umap and clusters to nucseq object - liger stores umap in tsne slot
  subset_so@reductions$umap@cell.embeddings <-
  	ligerex@tsne.coords[seurat.cell.ordering, ]
  dimnames(subset_so@reductions$umap@cell.embeddings)[[2]] <- c("UMAP_1", "UMAP_2")

  # Add liger clusters to seurat object
  subset_so <- AddMetaData(subset_so,
  	ligerex@alignment.clusters[seurat.cell.ordering],
  	col.name = "liger_clusters")

  return(subset_so)

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
###########################################################################
#
# find_cluster_enriched_genes_for_each_cluster <- function(
#   seurat_obj,
#   cluster_col_name){
#
#   print("find_cluster_enriched_genes_for_each_cluster")
#
#   clusters <- sort(unique(seurat_obj[[cluster_col_name]]))
#
#   top_expressed_genes_tb <-
#     map(clusters, function(cluster){
#       expr_m <- GetAssayData(seurat_obj, slot = "data")
#
#
#       # calculate percent of control cells expressing gene
#       dx_de_tb <- apply(expr_m, 1, function(row){
#           row_ctrl <- row[seurat_obj[["clinical_dx"]] == "Control"]
#           percent_detected <- (sum(row_ctrl > 0) / length(row_ctrl)) * 100
#           is.na(percent_detected) <- 0
#           return(percent_detected)
#         }) %>%
#         round(1) %>%
#         enframe(name = "gene", value = "percent_detected_ctrl") %>%
#         right_join(., dx_de_tb)
#
#       # calculate percent of dx cells expressing gene
#       dx_de_tb <- apply(expr_m, 1, function(row){
#           row_dx <- row[seurat_obj[["clinical_dx"]] != "Control"]
#           percent_detected <- (
#             sum(row_dx > 0) / length(row_dx)) * 100
#           is.na(percent_detected) <- 0
#           return(percent_detected)
#         }) %>%
#         round(1) %>%
#         enframe(name = "gene", value = "percent_detected_dx") %>%
#         right_join(., dx_de_tb)
#
#
#       [
#         , Idents(seurat_obj) == cluster] %>%
#       Matrix::rowMeans() %>%
#       sort(decreasing = TRUE) %>%
#       enframe() %>%
#       as_tibble() %>%
#       rename(gene = name, mean_expression = value) %>%
#       mutate(cluster = cluster)
#     }) %>%
#     bind_rows()
#
#   # subset seurat object to clinical_dxs and cell type of interest
#   # cell_type and clinical_dxs variables need to be global for some reason
#   seurat_obj <- subset(x = seurat_obj,
#     # subset = cluster_cell_type == cell_type_of_interest)
#     subset = cluster_cell_type_and_layer == cell_type_of_interest)
#   seurat_obj <- subset(x = seurat_obj,
#     subset = clinical_dx %in% clinical_dxs)
#   seurat_obj <- subset(x = seurat_obj,
#     subset = region == brain_region)
#
#   # DE Linear model
#   expr_m <- GetAssayData(seurat_obj, slot = "data")
#   terms_df <- data.frame(
#     seurat_obj[[c("clinical_dx", "number_umi"), drop = FALSE]])
#   terms_df$clinical_dx <- terms_df$clinical_dx %>%
#     factor %>%
#     relevel(., ref = "Control")
#   mod <- "y ~ clinical_dx + number_umi"
#   lm_coef_pval_l <- run_de_with_linear_model(
#     expr_m = expr_m, terms_df = terms_df, mod = mod)
#   dx_col_idx <- colnames(lm_coef_pval_l$coefmat) %>% grep("clinical_dx", .)
#
#   # format lm output into a tibble
#   dx_de_tb <- tibble(gene = rownames(lm_coef_pval_l$coefmat),
#     log2_fold_change_dx_vs_ctrl = lm_coef_pval_l$coefmat[ ,dx_col_idx],
#     pvalue = lm_coef_pval_l$pvalmat[ ,dx_col_idx]
#   )
#
#   # calculate percent of control cells expressing gene
#   dx_de_tb <- apply(expr_m, 1, function(row){
#       row_ctrl <- row[seurat_obj[["clinical_dx"]] == "Control"]
#       percent_detected <- (sum(row_ctrl > 0) / length(row_ctrl)) * 100
#       is.na(percent_detected) <- 0
#       return(percent_detected)
#     }) %>%
#     round(1) %>%
#     enframe(name = "gene", value = "percent_detected_ctrl") %>%
#     right_join(., dx_de_tb)
#
#   # calculate percent of dx cells expressing gene
#   dx_de_tb <- apply(expr_m, 1, function(row){
#       row_dx <- row[seurat_obj[["clinical_dx"]] != "Control"]
#       percent_detected <- (
#         sum(row_dx > 0) / length(row_dx)) * 100
#       is.na(percent_detected) <- 0
#       return(percent_detected)
#     }) %>%
#     round(1) %>%
#     enframe(name = "gene", value = "percent_detected_dx") %>%
#     right_join(., dx_de_tb)
#
#   # add ensembl
#   dx_de_tb$ensembl <- convert_gene_symbols_to_ensembl_ids(
#     gene_symbols = dx_de_tb$gene)
#
#   # fdr_pvalue correct
#   dx_de_tb$fdr_pvalue <- p.adjust(dx_de_tb$pvalue, method = "BH")
#   # check
#   table(dx_de_tb$pvalue < 0.05)
#   table(dx_de_tb$fdr_pvalue < 0.05)
#
#   # formating
#   dx_de_tb <- dx_de_tb %>%
#     # order columns
#     select(ensembl, gene, log2_fold_change_dx_vs_ctrl, pvalue, fdr_pvalue,
#       percent_detected_dx, percent_detected_ctrl) %>%
#     # order by log2 fold change
#     arrange(desc(log2_fold_change_dx_vs_ctrl)) %>%
#     # convert factors to characters
#     mutate_if(is.factor, as.character) %>%
#     as_tibble()
#
#   dx_de_tb <- dx_de_tb %>% mutate(
#     clinical_dxs = paste0(clinical_dxs, collapse = "_"),
#     brain_region = brain_region,
#     cell_type_of_interest = cell_type_of_interest)
#
#   return(dx_de_tb)
# }
#
# ## Function: DE Linear model
# # terms_df:
# # ExpCondition RIN.y     Seq.PC1      Seq.PC2
# #   CP         8.4       0.04792498   -0.090448567
# #   CP         8.1       0.53502697   -0.287629654
# #   CP         8.0       0.18824922   -0.155651102
# #   VZ         8.4       0.02529722   -0.100858264
# #   VZ         8.7       0.45139297   0.856908177
# #   VZ         9.1       0.27861748   -0.248868277
# # mod: "y~ExpCondition+RIN.y+Seq.PC1+Seq.PC2"
# run_de_with_linear_model <- function(expr_m, terms_df, mod) {
#   print("run_de_with_linear_model")
#   lm_fit <- apply(expr_m, 1
#     , function(y) {
#       mod <- as.formula(mod)
#       lm(mod, data = terms_df)})
#   coefmat <- matrix(NA, nrow = nrow(expr_m)
#     , ncol = length(coef(lm_fit[[1]])))
#   pvalmat <- matrix(NA, nrow = nrow(expr_m)
#     , ncol = length(summary(lm_fit[[1]])[[4]][ ,4]))
#   colnames(coefmat) <- names(coef(lm_fit[[1]]))
#   rownames(coefmat) <- rownames(expr_m)
#   colnames(pvalmat) <- names(summary(lm_fit[[1]])[[4]][ ,4])
#   rownames(pvalmat) <- rownames(expr_m)
#   for (i in 1:nrow(expr_m)) {
#     if (i%%100 == 0) {cat(".")}
#     coefmat[i, ] <- coef(lm_fit[[i]])
#     pvalmat[i, ] <- summary(lm_fit[[i]])[[4]][ ,4]
#   }
#   lm_coef_pval_l <- list(coefmat = coefmat, pvalmat = pvalmat)
#   return(lm_coef_pval_l)
# }
#
# convert_gene_symbols_to_ensembl_ids <- function(
#   gene_symbols,
#   in_biomart_gene_info = "../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"){
#   # Converts vector of gene symbols to ensembl IDs
#   # gene symbols must be character format
#   # Need to load:
#   biomart_df <- read.csv(in_biomart_gene_info, header = TRUE)
#   print("convert_gene_symbols_to_ensembl_ids")
#   idx <- match(gene_symbols, biomart_df$hgnc_symbol)
#   ens <- biomart_df$ensembl_gene_id[idx]
#   gene_symbols[! is.na(ens)] <- as.character(ens[! is.na(ens)])
#   return(ens)
# }
###########################################################################

### Testing Seurat FindMarkers functions

# Results:
# FindMarkers() with "negative binomial" uses the glm.nb function
# avg_logFC = log(mean(counts cell type) + 1) - log(mean(counts all cells) + 1)

#
# load("/u/scratch/d/dpolioud/liger_subcluster/20200218/p1-5_c1-5/percent_mt_8/reg_umi/k25_lambda3/liger_subcluster_bvFTD_Control_preCG_excitatory_test.rdat")

# test <- subset(subset_so, downsample = 50)
# fm <- FindMarkers(test, test.use = "negbinom", ident.1 = 0)
#
# dat_count <- data.frame(
#   expression = FetchData(object = test, slot = "counts", vars = c("DCC"))[,1],
#   cluster_id = ifelse(Idents(object = test) == "0", "cluster", "all")
#   # ,
#   # number_umi = FetchData(object = test, vars = c("number_umi"))[,1]
# )
# dat_count %>%
#   group_by(cluster_id) %>%
#   # summarize(mean_expression = mean(expression)) %>%
#   summarize(mean_expression = log(mean(expression)+1)) %>%
#   summarize(mean_expression[2] - mean_expression[1]) %>%
#   pull()
#
# head(fm,20)
#
# require(MASS)
# glm <- glm.nb(expression ~ cluster_id, data = dat_count)
# summary(glm)


# dat <- data.frame(
#   expression = FetchData(object = test, slot = "data", vars = c("DCC"))[,1],
#   cluster_id = ifelse(Idents(object = test) == "0", "cluster", "all")
#   # ,
#   # number_umi = FetchData(object = test, vars = c("number_umi"))[,1]
# )
# dat %>%
#   group_by(cluster_id) %>%
#   summarize(mean_expression = mean(expression)) %>%
#   summarize(mean_expression[2] - mean_expression[1]) %>%
#   pull()
#
# lm_fit <- lm(expression ~ cluster_id, data = dat)
# summary(lm_fit)
###########################################################################

### Run

main_function()
###########################################################################
