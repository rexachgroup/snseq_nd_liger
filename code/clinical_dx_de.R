# Damon Polioudakis
# 2019-09-24
# Run differential expression by clinical dx

# Must load modules:
#  module load R/3.6.0

# Sample qsub
#   qsub -N dx_de -t 1-42 -l h_data=64G,h_rt=12:00:00 -m n qsub_r_script.sh -p clinical_dx_de.R
################################################################################

rm(list = ls())
set.seed(27)

require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")

sessionInfo()

## Command args to input to loop through args and run DE as array job
cmd_args <- commandArgs(trailingOnly = TRUE)
print(paste0("cmd_args: ", cmd_args))

## Inputs
in_seurat <- paste0(
  "../analysis/seurat/20191011/percent_mt_8/reg_umi/"
  , "p1-5_c1-5_filtered.rdat")
marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")

## Variables
script_name <- "clinical_dx_de.R"
date <- format(Sys.Date(), "%Y%m%d")
graph_subtitle <- "P1-5 C1-5"

## Outputs
out_table <- paste0(
  "../analysis/clinical_dx_de/", date
  , "/percent_mt_8/tables/clinical_dx_de_")
out_graph <- paste0(
  "../analysis/clinical_dx_de/", date
  , "/percent_mt_8/graphs/clinical_dx_de_")
print(out_graph)
print(out_table)

# Make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
################################################################################

main_function <- function(){

  print("main_function")

  # browser()

  load(in_seurat)
  Idents(nd_so) <- nd_so[["RNA_snn_res.0.6"]]

  nd_so <- assign_cell_type_to_cell(seurat_obj = nd_so)
  nd_so <- assign_cell_type_cluster(
    seurat_obj = nd_so, cluster_col_name = "RNA_snn_res.0.6")

  # plot_dim_reduction_colored_by_module_score(seurat_obj = nd_so)
  # plot_dim_reduction_colored_by_individual_variable(
  #   seurat_obj = nd_so, var_name = "cell_type")
  # plot_dim_reduction_colored_by_cell_type_assignment(seurat_obj = nd_so)
  # plot_dim_reduction_colored_by_cluster_cell_type_assignment(seurat_obj = nd_so)
  # plot_cell_type_by_cluster_stacked_bar_plot(
  #   seurat_obj = nd_so, cluster_col_name = "RNA_snn_res.0.6")

  # make tibble of de function args to loop through with array job
  run_de_args_tb <- expand.grid(
      clinical_dx = list(
        c("PSP-S", "Control"),
        c("AD", "Control"),
        c("bvFTD", "Control")
      ),
      cell_type = c(
        "astrocyte1",
        "endo2",
        "excitatory3",
        "inhibitory4",
        "micro5",
        "oligo6",
        "OPC7"
      ),
      region = c(
        "preCG",
        "calcarine"
      )
    ) %>%
    as_tibble()
  print(paste0("Dimensions of run_de_args_tb: ",
    paste0(dim(run_de_args_tb), collapse = " ")))
  # need to assign in global environment for seurat subset() function for some reason
  assign("clinical_dxs",
    run_de_args_tb$clinical_dx[[as.numeric(cmd_args[1])]],
    envir = .GlobalEnv)
  assign("cell_type_of_interest",
    run_de_args_tb$cell_type[[as.numeric(cmd_args[1])]],
    envir = .GlobalEnv)
  assign("brain_region",
    run_de_args_tb$region[[as.numeric(cmd_args[1])]],
    envir = .GlobalEnv)

  # calculate DE with LM
  dx_de_tb <- run_de_with_lm_by_clinical_dx_and_cell_type(
    seurat_obj = nd_so,
    clinical_dxs = clinical_dxs,
    brain_region = brain_region,
    cell_type_of_interest = cell_type_of_interest)

  # save DE table to csv
  out_dx_de_csv <- paste0(
      out_table, paste0(clinical_dxs, collapse = "_"), "_",
      brain_region, "_",
      cell_type_of_interest, ".csv")
  print(paste0("out csv path: ", out_dx_de_csv))
  write_csv(dx_de_tb, path = out_dx_de_csv)

  print("end of main_function")
}
################################################################################

run_de_with_lm_by_clinical_dx_and_cell_type <- function(
  seurat_obj,
  clinical_dxs = c("PSP-S", "Control"),
  cell_type_of_interest = "excitatory3",
  brain_region = "calcarine"){

  print("run_de_with_lm_by_clinical_dx_and_cell_type")

  print(cell_type_of_interest)
  print(clinical_dxs)
  print(brain_region)

  # subset seurat object to clinical_dxs and cell type of interest
  # cell_type and clinical_dxs variables need to be global for some reason
  seurat_obj <- subset(x = seurat_obj,
    subset = cluster_cell_type == cell_type_of_interest)
  seurat_obj <- subset(x = seurat_obj,
    subset = clinical_dx == clinical_dxs)
  seurat_obj <- subset(x = seurat_obj,
    subset = region == brain_region)

  # DE Linear model
  expr_m <- GetAssayData(seurat_obj, slot = "data")
  terms_df <- data.frame(
    seurat_obj[[c("clinical_dx", "number_umi"), drop = FALSE]])
  terms_df$clinical_dx <- terms_df$clinical_dx %>%
    factor %>%
    relevel(., ref = "Control")
  mod <- "y ~ clinical_dx + number_umi"
  lm_coef_pval_l <- run_de_with_linear_model(
    expr_m = expr_m, terms_df = terms_df, mod = mod)
  dx_col_idx <- colnames(lm_coef_pval_l$coefmat) %>% grep("clinical_dx", .)

  # format lm output into a tibble
  dx_de_tb <- tibble(gene = rownames(lm_coef_pval_l$coefmat),
    log2_fold_change_dx_vs_ctrl = lm_coef_pval_l$coefmat[ ,dx_col_idx],
    pvalue = lm_coef_pval_l$pvalmat[ ,dx_col_idx]
  )

  # calculate percent of control cells expressing gene
  dx_de_tb <- apply(expr_m, 1, function(row){
      row_ctrl <- row[seurat_obj[["clinical_dx"]] == "Control"]
      percent_detected <- (sum(row_ctrl > 0) / length(row_ctrl)) * 100
      is.na(percent_detected) <- 0
      return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_ctrl") %>%
    right_join(., dx_de_tb)

  # calculate percent of dx cells expressing gene
  dx_de_tb <- apply(expr_m, 1, function(row){
      row_dx <- row[seurat_obj[["clinical_dx"]] != "Control"]
      percent_detected <- (
        sum(row_dx > 0) / length(row_dx)) * 100
      is.na(percent_detected) <- 0
      return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_dx") %>%
    right_join(., dx_de_tb)

  # add ensembl
  dx_de_tb$ensembl <- convert_gene_symbols_to_ensembl_ids(
    gene_symbols = dx_de_tb$gene)

  # fdr_pvalue correct
  dx_de_tb$fdr_pvalue <- p.adjust(dx_de_tb$pvalue, method = "BH")
  # check
  table(dx_de_tb$pvalue < 0.05)
  table(dx_de_tb$fdr_pvalue < 0.05)

  # formating
  dx_de_tb <- dx_de_tb %>%
    # order columns
    select(ensembl, gene, log2_fold_change_dx_vs_ctrl, pvalue, fdr_pvalue,
      percent_detected_dx, percent_detected_ctrl) %>%
    # order by log2 fold change
    arrange(desc(log2_fold_change_dx_vs_ctrl)) %>%
    # convert factors to characters
    mutate_if(is.factor, as.character) %>%
    as_tibble()

  dx_de_tb <- dx_de_tb %>% mutate(
    clinical_dxs = paste0(clinical_dxs, collapse = "_"),
    brain_region = brain_region,
    cell_type_of_interest = cell_type_of_interest)

  return(dx_de_tb)
}

## Function: DE Linear model
# terms_df:
# ExpCondition RIN.y     Seq.PC1      Seq.PC2
#   CP         8.4       0.04792498   -0.090448567
#   CP         8.1       0.53502697   -0.287629654
#   CP         8.0       0.18824922   -0.155651102
#   VZ         8.4       0.02529722   -0.100858264
#   VZ         8.7       0.45139297   0.856908177
#   VZ         9.1       0.27861748   -0.248868277
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

convert_gene_symbols_to_ensembl_ids <- function(
  gene_symbols,
  in_biomart_gene_info = "../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"){
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
################################################################################


main_function()
################################################################################
