set.seed(27)
require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
require(future.apply)
require(ggpubr)
require(patchwork)
require(lme4)
require(lmerTest)
require(broom.mixed)
require(ComplexHeatmap)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
options(future.globals.maxSize = Inf)

# in_seurat_expr <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_expression_matrix.rdat"
# in_seurat_meta <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_metadata.csv"
in_seurat_rds <-
  "../analysis/pci_import/pci_seurat.rds"

marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
marker_genes_refined_tb <- read_csv(
  "../resources/20191028_cell_markers_refined.csv")
in_exc_marker <- "../resources/excitatory_markers_20191023.csv"
in_jessica_go_gene_list <- "../resources/20200607_jessica_input_up.txt"
metadata_corrections_tb <- read_csv("../metadata/2020.06.13_corrected_metadata_jessica.csv")

## Variables
script_name <- "clinical_dx_de_microglia.R"
date <- format(Sys.Date(), "%Y%m%d")
graph_subtitle <- "P1-5"

## Outputs
out_path_base <- "../analysis/clinical_dx_de_lme_microglia/"
out_table <- paste0(
  out_path_base, date,
  "/tables/clinical_dx_de_lme_microglia_")
out_graph <- paste0(
  out_path_base, date,
  "/graphs/clinical_dx_de_lme_microglia_")
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_graph), recursive = TRUE)
print(out_graph)
print(out_table)

main <- function() {
    nd_so <- readRDS(in_seurat_rds)

    #     clinical_dx <- c("bvFTD", "Control")
    cell_type <- "microglia"
    region <- "preCG"

    # Filter metadata to get cell ids for matching region, dx, cluster_cell_type.
    # Subset seurat object.
    cell_id_subset <- nd_so@meta.data %>%
        filter(cell_type == {{cell_type}},
               region == {{region}}) %>%
        pluck("cell_ids")
    mc_so <- subset(nd_so, cells = cell_id_subset)
    mc_so <- nd_so[, cell_id_subset]
    meta <- mc_so@meta.data
    rm(nd_so)
    gc()
    mc_so[["log_number_umi"]] <- log(mc_so[["number_umi"]])
    
    # Run model.
    model_design <- "expression ~ clinical_dx + pmi + age + rin + sex + seq_batch + number_umi + percent_mito + Reads.Mapped.Antisense.to.Gene + Fraction.Reads.in.Cells + (1 | library_id)"
    plan(multicore, workers = 12)
    lm_broom <- run_lmer_de(mc_so, model_design, down_sample_cells = 10000)
    print(sum(is.na(lm_broom)))
    plan(multicore)
    lm_tb <- filter_and_format_lm_output(lm_broom[!is.na(lm_broom)], mc_so, model_design, beta_regex = "clinical_dx")
    plan(sequential)


    # Write out lm_broom. Format and write out lm_tb.
    saveRDS(lm_broom, paste0(dirname(out_table), "/lm.rds"))
    saveRDS(lm_tb, paste0(out_table, "lm_tb.rds"))
    write_csv(lm_tb, paste0(out_table, "lm_tb.csv"))
    
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = -c("gene", "model", "percent_detected_ctrl", "percent_detected_dx", "term"),
        names_glue = "{term}.{.value}"
    ) %>% arrange(term)

    lm_tb %>% 
        pivot_wider_spec(lm_wider_spec) %>%
        write_csv(., paste0(out_table, "lm_wider_tb.csv"))    
}


run_lmer_de <- function(
    seurat_obj,
    model_design,
    down_sample_cells = NULL,
    cores = NULL){
 
    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    # If dx is present relevel seurat_obj that control is the 1st factor level.
    model <- as.formula(model_design)
    test_vars <- seurat_obj@meta.data

    if ("clinical_dx" %in% colnames(test_vars)) {
        test_vars <- mutate(test_vars, clinical_dx = fct_relevel(clinical_dx, "Control"))
    }
    
    run_lmer_broom(expr_m, test_vars, model)
}

run_lmer_broom <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr_r) {
        dat <- data.frame(expression = as.vector(expr_r), test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lmer(model, data = dat))
        },
        error = function(x) { print(x); NA })
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

    if (nrow(lm_out_tb) == 0)
        stop("lm output is of length zero. Check lm_out_obj_l and beta_regex arguments")

    expr_m <- GetAssayData(seurat_obj, slot = "data")
    is_ctrl <- seurat_obj@meta.data[["clinical_dx"]] == "Control"
    is_dx <- seurat_obj@meta.data[["clinical_dx"]] != "Control"
    expr_in_lmtb <- expr_m[unique(lm_out_tb$gene), ]

    # percent of control cells expressing gene
    percent_detected_ctrl <- future_apply(expr_in_lmtb, 1, function(row){
        row_ctrl <- row[is_ctrl]
        percent_detected <- (sum(row_ctrl > 0) / length(row_ctrl)) * 100
        percent_detected[!is.finite(percent_detected)] <- 0
        return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_ctrl")

    # percent of dx cells expressing gene
    percent_detected_dx <- future_apply(expr_in_lmtb, 1, function(row){
        row_dx <- row[is_dx]
        percent_detected <- (sum(row_dx > 0) / length(row_dx)) * 100
        percent_detected[!is.finite(percent_detected)] <- 0
        return(percent_detected)
    }) %>%
    round(1) %>%
    enframe(name = "gene", value = "percent_detected_dx")

    lm_out_tb <- left_join(lm_out_tb, percent_detected_ctrl, by = "gene") %>%
        left_join(percent_detected_dx, by = "gene")

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
    arrange(desc(estimate))
    # convert factors to characters
  # add additional columns
  lm_out_tb[["model"]] <- model_design

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

ggplot_heatmap_cluster <- function(tb_cor) {
    
}

if (interactive()) {
    main()
}
