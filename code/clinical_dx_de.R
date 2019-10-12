# Damon Polioudakis
# 2019-09-24
# Run differential expression by clinical dx

# Must load modules:
#  module load R/3.6.0

# Sample qsub
#   qsub -N dx_de -t 1-42 -l h_data=64G,h_rt=12:00:00 qsub_r_script.sh -p clinical_dx_de.R
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
cmd_args <- 38
print(paste0("cmd_args: ", cmd_args))

## Inputs
in_seurat <- paste0(
  "../analysis/seurat/20190916/percent_mt_8/reg_mt/"
  , "p1-5_c1-5_filtered_rmp27.rdat")
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

assign_cell_type_to_cell <- function(seurat_obj){

  print("assign_cell_type_to_cell")

  marker_genes <- marker_genes_tb %>%
    filter(source == "tsai") %>%
    filter(gene_symbol %in% rownames(seurat_obj)) %>%
    split(x = ., f = .$marker_for) %>%
    map(., "gene_symbol")

  # AddModuleScore() appends a sequential digit to the end of the label
  # e.g. astrocyte to astrocyte1, endo to endo2
  seurat_obj <- AddModuleScore(object = seurat_obj, features = marker_genes,
    name = names(marker_genes))

  seurat_obj$cell_type <- seurat_obj[[c("cell_ids", "astrocyte1", "endo2", "excitatory3", "inhibitory4", "micro5", "oligo6", "OPC7")]] %>%
    as_tibble() %>%
    gather(key = "cell_type", value = "score", -cell_ids) %>%
    group_by(cell_ids) %>%
    slice(which.max(score)) %>%
    select(cell_ids, cell_type) %>%
    right_join(., seurat_obj[["cell_ids"]]) %>%
    pull(cell_type)

  return(seurat_obj)

}

assign_cell_type_cluster <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6"){

  print("assign_cell_type_cluster")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]] %>% pull

  seurat_obj$cluster_cell_type <-
    seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]] %>%
      as_tibble() %>%
      group_by(cluster_ids) %>%
      count(cell_type) %>%
      mutate(percent = n/sum(n)*100) %>%
      slice(which.max(percent)) %>%
      mutate(cluster_cell_type = ifelse(percent > 50, cell_type, "mixed")) %>%
      right_join(., seurat_obj[[c("cluster_ids", "cell_ids", "cell_type")]],
        by = "cluster_ids") %>%
      pull(cluster_cell_type)

  return(seurat_obj)

}

plot_dim_reduction_colored_by_module_score <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_module_score")

  cell_marker_scores_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "astrocyte1", "endo2", "excitatory3", "inhibitory4", "micro5", "oligo6", "OPC7")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

    gg_l <- map(c("RNA_snn_res.0.6", "astrocyte1", "endo2", "excitatory3", "inhibitory4", "micro5", "oligo6", "OPC7"), function(variable){
      gg_tb <- cell_marker_scores_tb
      plot_dim_reduction_colored_by_variable(
        dim_1 = gg_tb$dim_1
        , dim_2 = gg_tb$dim_2
        , variable_value = gg_tb[[variable]]
        , title = variable
        , legend_title = variable
        , expression_color_gradient = TRUE
        , size = 0.5
        , alpha = 0.5
        , guide_size = 2
      )
    })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.3
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by marker expression score", "\n", reduction))
    )
  ggsave(paste0(out_graph, "marker_score_", reduction, ".png")
    , width = 20, height = 8, dpi = 200)

}

plot_dim_reduction_colored_by_individual_variable <- function(
  seurat_obj, var_name = "cell_type", reduction = "umap",
  factor_to_numeric = FALSE){

  print("plot_dim_reduction_colored_by_individual_variable")

  seurat_obj$variable_id <- seurat_obj[[var_name]]

  var_dim_reduction_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[["variable_id"]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble() %>%
      { if(factor_to_numeric == TRUE){
        mutate(., variable_id = variable_id %>%
          as.character() %>%
          as.numeric())
      } else .
      }
    )
  variables <- var_dim_reduction_tb$variable_id %>%
    unique() %>% sort()
  gg_l <- map(variables, function(variable){
    var_dim_reduction_tb %>%
      mutate(variable_membership = if_else(
        variable_id == variable, TRUE, FALSE)) %>%
      # sort so that TRUE points are printed on top of FALSE
      arrange(variable_membership) %>%
        # pull(variable_membership) %>% table
      {plot_dim_reduction_colored_by_variable(
        dim_1 = .$dim_1
        , dim_2 = .$dim_2
        , variable_value = .$variable_membership
        , title = paste0(var_name, ": ", variable)
        , legend_title = paste0(var_name, ": ", variable)
        , size = 0.1
        , alpha = 0.5
        , guide_size = 2
      ) + scale_color_manual(values = c("grey", "#00BFC4"))}
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.3
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by ", var_name, "\n", reduction))
    )
  ggsave(paste0(out_graph, "variable_facet_", var_name, "_", reduction, ".png")
    , width = 18, height = 3.25+0.6*length(variables), dpi = 200, limitsize = FALSE)
}

plot_cell_type_by_cluster_stacked_bar_plot <- function(
  seurat_obj, cluster_col_name = "RNA_snn_res.0.6"){

  print("plot_cell_type_by_cluster_stacked_bar_plot")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # collect metadata
  metadata_tb <- seurat_obj[[]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

    # percent stacked bar chart
    # calculate percent
    gg_tb <- metadata_tb %>%
      filter(variable == "cell_type") %>%
      group_by(cluster_ids) %>%
      count(value) %>%
      mutate(percent = n/sum(n)) %>%
      mutate(variable = "cell_type") %>%
      ungroup() %>%
      mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
      mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))
    # plot
    ggplot(gg_tb, aes(x = cluster_ids, y = percent, fill = value)) +
      geom_bar(stat = "identity") +
      { if (length(unique(gg_tb$value)) > 11){
        scale_fill_manual(name = "cell_type", values = colorRampPalette(
          brewer.pal(n = 11, name = "Set3"))(length(unique(gg_tb$value))))
      } else {
        scale_fill_brewer(name = "cell_type", palette = "Set3")
      }} +
      ggtitle(make_plot_title("Percent assigned cell type by cluster"))
    ggsave(paste0(out_graph, "cell_type_by_cluster_percent_stackedbar.png"),
      width = 12, height = 6, limitsize = FALSE)

}

plot_dim_reduction_colored_by_cell_type_assignment <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_cell_type_assignment")

  cell_type_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "cell_type")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

  plot_dim_reduction_colored_by_variable(
    dim_1 = cell_type_tb$dim_1
    , dim_2 = cell_type_tb$dim_2
    , variable_value = cell_type_tb[["cell_type"]]
    , title = make_plot_title(paste0("Cell type assigned by max cell type marker score", "\n", reduction))
    , legend_title = "cell_type"
    , size = 0.5
    , alpha = 0.5
    , guide_size = 2
  )

  ggsave(paste0(out_graph, "cell_type_", reduction, ".png")
    , width = 7, height = 6, dpi = 200)

}

plot_dim_reduction_colored_by_cluster_cell_type_assignment <- function(
  seurat_obj, reduction = "umap"){

  print("plot_dim_reduction_colored_by_cluster_cell_type_assignment")

  cell_type_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble() %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[[c("RNA_snn_res.0.6", "cluster_cell_type")]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble()
    )

  plot_dim_reduction_colored_by_variable(
    dim_1 = cell_type_tb$dim_1
    , dim_2 = cell_type_tb$dim_2
    , variable_value = cell_type_tb[["cluster_cell_type"]]
    , title = make_plot_title(paste0("Cell type assigned to cluster by cell type marker score", "\n", reduction))
    , legend_title = "cluster_cell_type"
    , size = 0.5
    , alpha = 0.5
    , guide_size = 2
  )

  ggsave(paste0(out_graph, "cluster_cell_type_", reduction, ".png")
    , width = 7, height = 6, dpi = 200)

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
    cell_type_of_interest = cell_type_of_interest,
    number_of_cells = number_of_cells)

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
