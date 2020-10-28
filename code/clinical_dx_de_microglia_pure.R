set.seed(27)
require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
require(broom)
require(future.apply)
require(ggpubr)
require(patchwork)
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
script_name <- "clinical_dx_de_microglia_pure.R"
date <- format(Sys.Date(), "%Y%m%d")
graph_subtitle <- "P1-5"

## Outputs
out_path_base <- "../analysis/clinical_dx_de_microglia_pure/"
out_table <- paste0(
  out_path_base, date,
  "/tables/clinical_dx_de_microglia_pure_")
out_graph <- paste0(
  out_path_base, date,
  "/graphs/clinical_dx_de_microglia_pure_")
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_graph), recursive = TRUE)
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

    # Filter metadata to get cell ids for matching region, dx, cluster_cell_type.
    # Subset seurat object.
    cell_id_subset <- so@meta.data %>%
        filter(region == {{region}},
               clinical_dx %in% {{clinical_dx}},
               cell_type_and_layer == {{cell_type}}) %>%
        rownames
    mc_so <- subset(so, cells = cell_id_subset)

    
    model_design <- "expression ~ clinical_dx + sex + finalsite + pmi_h + age + rin + percent_mito + prep"
    lm_broom <- run_lm_de(
        seurat_obj = mc_so,
        model_design = model_design,
        down_sample_cells = 10000)
    lm_tb <- filter_and_format_lm_output(
        lm_out_obj_l = lm_broom, 
        seurat_obj = mc_so,
        model_design = model_design,
        beta_regex = "clinical_dx")

    precg_micro_bvftd <- mc_so@meta.data %>%
        as_tibble(rownames = "id")

    test_tb <- wilcox_test_cols(precg_micro_bvftd, "clinical_dx", 
        list("age", "rin", "percent_mito", "pmi_h", "number_umi", "number_genes", "finalsite", "sex", "prep"))

    meta_barplots <- ggplot_grouped_boxplot(precg_micro_bvftd, "clinical_dx",
        c("age", "rin", "percent_mito", "pmi_h", "number_umi", "number_genes", "prep"))

    write_csv(precg_micro_bvftd, paste0(out_table, "meta.csv"))
    write_csv(lm_tb, paste0(out_table, "lm_tb.csv"))
    write_csv(test_tb, paste0(out_table, "mann-whitney_metadata.csv"))

    pdf(paste0(out_graph, "meta_barplots.pdf"), width = 16, height = 8 * (length(meta_barplots) / 2))
    tryCatch({wrap_plots(meta_barplots, ncol = 2)}, error = print)
    dev.off()
}


run_lm_de <- function(
    seurat_obj,
    model_design,
    down_sample_cells = NULL,
    cores = NULL){
 
    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        seurat_obj <- subset(seurat_obj, cells = sample(rownames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    # Subset metadata to formula terms.
    # If dx is present relevel seurat_obj that control is the 1st factor level.
    model <- as.formula(model_design)
    test_vars <- seurat_obj@meta.data %>%
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
    model_design,
    beta_regex,
    percent_detected_filter = 0){

  print("filter_and_format_lm_output")

  clinical_dx <- unique(seurat_obj[["clinical_dx"]])
  region <- unique(seurat_obj[["region"]])
  cell_type <- unique(seurat_obj[["cell_type"]])
  # format lm output into tibble
  lm_out_tb <- bind_rows(lm_out_obj_l, .id = "gene") %>%
    filter(grepl(beta_regex, term) | is.na(term))

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
#   lm_out_tb$clinical_dx <- paste0(clinical_dx, collapse = "_")
#   lm_out_tb$region <- region
#   lm_out_tb$cell_type <- cell_type
#   lm_out_tb$model_design <- model_design

  print("end of filter_and_format_lm_output")

  return(lm_out_tb)

}


wilcox_test_cols <- function(tb, group_name, test_list) {
    lapply(test_list, function(var_name) {
        test_vec <- tb[[var_name]]
        if (!is.numeric(test_vec)) {
            test_vec <- as.numeric(factor(test_vec))
            var_name <- paste0("as_numeric_", var_name)
        }
        group_vec <- tb[[group_name]]
        out_tb <- broom::tidy(pairwise.wilcox.test(test_vec, group_vec, paired = FALSE))
        out_tb$var <- var_name
        return(out_tb)
    }) %>% bind_rows()
}

ggplot_grouped_boxplot <- function(tb, x_name, var_list) {
    return(lapply(var_list, function(var_name) {
        ggplot(tb, aes_string(x = x_name, y = var_name, fill = x_name)) +
            geom_boxplot() +
            stat_compare_means(paired = FALSE) +
            ggtitle(var_name)
    }))
}

if (interactive()) {
    main()
}
