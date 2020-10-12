set.seed(27)
require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
require(broom)
require(future.apply)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
plan(multicore)

in_seurat_expr <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_expression_matrix.rdat"
in_seurat_meta <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_metadata.csv"

marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
marker_genes_refined_tb <- read_csv(
  "../resources/20191028_cell_markers_refined.csv")
in_exc_marker <- "../resources/excitatory_markers_20191023.csv"
in_jessica_go_gene_list <- "../resources/20200607_jessica_input_up.txt"
metadata_corrections_tb <- read_csv("../metadata/2020.06.13_corrected_metadata_jessica.csv")

## Variables
script_name <- "clinical_dx_de_lm_pmi_brain_bank.R"
date <- format(Sys.Date(), "%Y%m%d")
graph_subtitle <- "P1-5"

## Outputs
out_path_base <- "../analysis/clinical_dx_de_pmi_brain_bank/"
out_table <- paste0(
  out_path_base, date,
  "/tables/clinical_dx_de_pmi_brain_bank_")
out_graph <- paste0(
  out_path_base, date,
  "/graphs/clinical_dx_de_pmi_brain_bank_")
seurat_checkpt <- paste0(out_path_base, date, "/p1-5_c1-5_fixed.rds")
print(out_graph)
print(out_table)

main <- function() {
    load(in_seurat_expr)
    meta <- read_csv(in_seurat_meta) %>%
        as.data.frame() %>%
        column_to_rownames("cell_id")

    so <- CreateSeuratObject(counts = normalized_expression_matrix, meta.data = meta)
    model_design_dx <- "expression ~ clinical_dx + number_umi"
    so <- assign_excitatory_layer_to_cell(seurat_obj = so)
    so <- assign_cell_type_and_layer_to_cluster(
      seurat_obj = so, cluster_col_name = "cluster_ids")


    clinical_dx <- c("bvFTD", "Control")
    cell_type <- "microglia"
    region <- "preCG"
    # calculate DE with LM
    model_design_set1 <- "expression ~ clinical_dx + sex + finalsite + pmi_h"
    lm_broom_set1 <- run_de_with_lm_by_clinical_dx_and_cell_type(
        seurat_obj = so,
        clinical_dx = clinical_dx,
        region = region,
        cell_type = cell_type,
        model_design = model_design_set1,
        down_sample_cells = 10000)
    lm_tb_set1 <- filter_and_format_lm_output(
        lm_out_obj_l = lm_broom_set1, 
        seurat_obj = so,
        clinical_dx = clinical_dx,
        region = region,
        cell_type = cell_type,
        model_design = model_design_set1,
        beta_regex = "clinical_dx")
    
    model_design_set2 <- "expression ~ clinical_dx + sex + finalsite + pmi_h + age + rin + percent_mito"
    lm_broom_set2 <- run_de_with_lm_by_clinical_dx_and_cell_type(
        seurat_obj = so,
        clinical_dx = clinical_dx,
        region = region,
        cell_type = cell_type,
        model_design = model_design_set2,
        down_sample_cells = 10000)
    lm_tb_set2 <- filter_and_format_lm_output(
        lm_out_obj_l = lm_broom_set2, 
        seurat_obj = so,
        clinical_dx = clinical_dx,
        region = region,
        cell_type = cell_type,
        model_design = model_design_set2,
        beta_regex = "clinical_dx")


    model_pmi <- "expression ~ pmi_h"
    lm_broom_pmi_ctl <- run_de_with_lm_by_clinical_dx_and_cell_type(
        seurat_obj = so,
        clinical_dx = "Control",
        region = region,
        cell_type = cell_type,
        model_design = model_pmi,
        down_sample_cells = 10000)
    lm_broom_pmi_dx_bvftd <- run_de_with_lm_by_clinical_dx_and_cell_type(
        seurat_obj = so,
        clinical_dx = "bvFTD",
        region = region,
        cell_type = cell_type,
        model_design = model_pmi,
        down_sample_cells = 10000)

    lm_broom_pmi_dx_pid <- run_de_with_lm_by_clinical_dx_and_cell_type(so, "PiD", region, cell_type, model_pmi, 10000)

    lm_tb_pmi_ctl <- filter_and_format_lm_output(lm_broom_pmi_ctl, so, 
        clinical_dx = "Control", cell_type = cell_type, region = region, model = model_pmi, "pmi_h")
    lm_tb_pmi_dx_bvftd <- filter_and_format_lm_output(lm_broom_pmi_dx_bvftd, so, 
        clinical_dx = "bvFTD", cell_type = cell_type, region = region, model_pmi, "pmi_h")


    write_csv(lm_tb_set1, paste0(out_table, "lm_tb_set1.csv"))
    write_csv(lm_tb_set2, paste0(out_table, "lm_tb_set2.csv"))
    write_csv(lm_tb_pmi_ctl, paste0(out_table, "lm_tb_pmi_ctl.csv"))
    write_csv(lm_tb_pmi_dx, paste0(out_table, "lm_tb_pmi_dx.csv"))
}


run_de_with_lm_by_clinical_dx_and_cell_type <- function(
    seurat_obj,
    clinical_dx,
    cell_type,
    region,
    model_design,
    down_sample_cells = NULL,
    cores = NULL){

    # Filter metadata to get cell ids for matching region, dx, cluster_cell_type.
    # Subset seurat object.
    cell_id_subset <- seurat_obj@meta.data %>%
        filter(region == {{region}},
               clinical_dx %in% {{clinical_dx}},
               cluster_cell_type_and_layer == {{cell_type}}) %>%
        rownames

    so <- subset(seurat_obj, cells = cell_id_subset)
    
    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(so) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        so <- subset(so, cells = sample(rownames(so), down_sample_cells))
        so@meta.data %>%
            select(region, clincal_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(so, slot = "data")

    # Convert model to formula.
    # Subset metadata to formula terms.
    # If dx is present relevel so that control is the 1st factor level.
    model <- as.formula(model_design)
    test_vars <- so@meta.data %>%
        select(all_of(attr(terms(model), "term.labels")))

    if ("clinical_dx" %in% colnames(test_vars)) {
        test_vars <- mutate(test_vars, clinical_dx = fct_relevel(clinical_dx, "Control"))
    }
    
    run_lm_broom(expr_m, test_vars, model)
}

run_lm_broom <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr) {
        dat <- data.frame(expression = expr, test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lm(model, data = dat))
        },
        error = function(e){NA})
    })
    return(lm_out_obj_l)
}

filter_and_format_lm_output <- function(
    lm_out_obj_l, 
    seurat_obj,
    clinical_dx,
    cell_type,
    region,
    model_design,
    beta_regex,
    percent_detected_filter = 0){

  print("filter_and_format_lm_output")

  # format lm output into tibble
  lm_out_tb <- bind_rows(lm_out_obj_l, .id = "gene") %>%
    filter(grepl(beta_regex, term) | is.na(term))

  # calculate percent of cells expressing gene
  # subset seurat object to clinical_dx and cell type of interest
    cell_id_subset <- seurat_obj@meta.data %>%
        filter(region == {{region}},
               clinical_dx %in% {{clinical_dx}},
               cluster_cell_type_and_layer == {{cell_type}}) %>%
        rownames
    seurat_obj <- subset(seurat_obj, cells = cell_id_subset)
  expr_m <- GetAssayData(seurat_obj, slot = "data")
  # percent of control cells expressing gene
  lm_out_tb <- apply(expr_m, 1, function(row){
      row_ctrl <- row[seurat_obj[["clinical_dx"]] == "Control"]
      percent_detected <- (sum(row_ctrl > 0) / length(row_ctrl)) * 100
      is.na(percent_detected) <- 0
      return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_ctrl") %>%
    right_join(., lm_out_tb)
  # percent of dx cells expressing gene
  lm_out_tb <- apply(expr_m, 1, function(row){
      row_dx <- row[seurat_obj[["clinical_dx"]] != "Control"]
      percent_detected <- (
        sum(row_dx > 0) / length(row_dx)) * 100
      is.na(percent_detected) <- 0
      return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_dx") %>%
    right_join(., lm_out_tb)

  # filter to genes expressed in greater than X percent of cells
  print(paste0("number of genes before percent detected filter: ", nrow(lm_out_tb)))
  lm_out_tb <- lm_out_tb %>%
    filter(percent_detected_dx > percent_detected_filter | percent_detected_ctrl > percent_detected_filter)
  print(paste0("number of genes after percent detected filter: ", nrow(lm_out_tb)))

  # add ensembl
  #   lm_out_tb$ensembl <- convert_gene_symbols_to_ensembl_ids(
  #     gene_symbols = lm_out_tb$gene)

  # fdr_pvalue correct
  lm_out_tb$fdr_pvalue <- p.adjust(lm_out_tb$p.value, method = "BH")
  # check
  table(lm_out_tb$p.value < 0.05) %>% print()
  table(lm_out_tb$fdr_pvalue < 0.05) %>% print()

  # formating
  lm_out_tb <- lm_out_tb %>%
    # order columns
    dplyr::select(
      gene,
      term = term,
      estimate = estimate,
      t_statistic = statistic,
      percent_detected_dx,
      percent_detected_ctrl,
      pvalue = p.value,
      fdr_pvalue,
      everything()) %>%
    clean_variable_names() %>%
    # order by log2 fold change
    arrange(desc(estimate)) %>%
    # convert factors to characters
    mutate_if(is.factor, as.character) %>%
    as_tibble()
  # add additional columns
  lm_out_tb <- lm_out_tb %>% mutate(
    clinical_dx = paste0(clinical_dx, collapse = "_"),
    region = region,
    cell_type = cell_type,
    model_design = model_design)

  print("end of filter_and_format_lm_output")

  return(lm_out_tb)

}

if (interactive()) {
    main()
}
