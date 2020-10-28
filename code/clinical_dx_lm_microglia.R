set.seed(27)
require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(RColorBrewer)
require(broom)
require(broom.mixed)
require(future.apply)
require(ggpubr)
require(ggfortify)
require(patchwork)
require(nestedRanksTest)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
options(future.globals.maxSize = Inf)

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
script_name <- "clinical_dx_de_microglia.R"
date <- format(Sys.Date(), "%Y%m%d")
graph_subtitle <- "P1-5"

## Outputs
out_path_base <- "../analysis/clinical_dx_lm_microglia/"
out_table <- paste0(
  out_path_base, date,
  "/tables/clinical_dx_lm_microglia")
out_graph <- paste0(
  out_path_base, date,
  "/graphs/clinical_dx_lm_microglia")
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_graph), recursive = TRUE)
print(out_graph)
print(out_table)

main <- function() {
    plan(multicore)
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
    
    cell_id_subset <- so@meta.data %>%
        filter(region == {{region}},
               clinical_dx %in% {{clinical_dx}},
               cluster_cell_type_and_layer == {{cell_type}}) %>%
        rownames
    mc_so <- subset(so, cells = cell_id_subset)

    # dx
    model_dx <- "expression ~ clinical_dx"
    lm_list_dx <- run_lm_de(mc_so, model_dx, down_sample_cells = 10000)
    lm_tb_dx <- filter_and_format_lm_output(lm_list_dx, mc_so, model_dx, "clinical_dx")
    
    # dx + umi.
    model_dx_umi <- "expression ~ clinical_dx + number_umi"
    lm_list_dx_umi <- run_lm_de(mc_so, model_dx_umi, down_sample_cells = 10000)
    lm_tb_dx_umi <- filter_and_format_lm_output(lm_list_dx_umi, mc_so, model_dx_umi, "clinical_dx")

    # dx, umi, percent_mito.
    model_dx_umi_mito <- "expression ~ clinical_dx + number_umi + percent_mito"
    lm_list_dx_umi_mito <- run_lm_de(mc_so, model_dx_umi_mito, down_sample_cells = 10000)
    lm_tb_dx_umi_mito <- filter_and_format_lm_output(lm_list_dx_umi_mito, mc_so, model_dx_umi_mito, "clinical_dx")
    
    # Full model.
    model_design <- "expression ~ clinical_dx + sex + finalsite + pmi_h + age + rin + percent_mito + prep + number_umi"
    lm_list <- run_lm_de(mc_so, model_design, down_sample_cells = 10000)
    lm_tb <- filter_and_format_lm_output(lm_list, mc_so, model_design, beta_regex = "clinical_dx") 
    
    # per-variable models.
    test_vars <- c("sex", "finalsite", "pmi_h", "age", "rin", "prep", "number_genes")
    individual_vars <- map(test_vars, function(test_var) {
        map(clinical_dx, function(dx) {
            writeLines(str_glue("========== {test_var}, {dx} =========="))
            cell_id_subset <- mc_so@meta.data %>%
                filter(clinical_dx == {{dx}}) %>%
                rownames
            test_so <- subset(mc_so, cells = cell_id_subset)
            model <- str_glue("expression ~ {test_var} + percent_mito + number_umi")
            i_lm_list <- run_lm_de(test_so, model, down_sample_cells = 10000)
            i_lm_tb <- filter_and_format_lm_output(i_lm_list, test_so, model, beta_regex = str_glue("^{test_var}"))

            i_lm_tb$model <- model
            i_lm_tb$dx <- dx

            i_lm_tb <- select(i_lm_tb, model, dx, everything())
            return(i_lm_tb)
        })
    })

    # partial models.
    lm_pmi_age_site_sex_design <- "expression ~ clinical_dx + pmi_h + age + finalsite + sex"
    lm_pmi_age_site_sex_list <- run_lm_de(mc_so, lm_pmi_age_site_sex_design, down_sample_cells = 10000)
    lm_pmi_age_site_sex_tb <- filter_and_format_lm_output(lm_pmi_age_site_sex_list, mc_so, lm_pmi_age_site_sex_design, beta_regex = "clinical_dx")
    
    lm_pmi_pct_mito_design <- "expression ~ clinical_dx + sex + finalsite + pmi_h + age + rin + number_umi + number_genes"
    lm_pmi_pct_mito_list <- run_lm_de(mc_so, lm_pmi_pct_mito_design, down_sample_cells = 10000)
    lm_pmi_pct_mito_tb <- filter_and_format_lm_output(lm_pmi_pct_mito_list, mc_so, lm_pmi_pct_mito_design, beta_regex = "clinical_dx")

    lm_ctl_rin_umi_design <- "expression ~ rin + number_umi"
    ctl_cells <- mc_so@meta.data %>% filter(clinical_dx == "Control") %>% rownames
    ctl_so <- subset(mc_so, cells = ctl_cells)
    lm_ctl_rin_umi_tb <- run_lm_de(ctl_so, lm_ctl_rin_umi_design, 10000) %>%
        filter_and_format_lm_output(ctl_so, lm_ctl_rin_umi_design, beta_regex = "^[^(]")

    # wilcox variable testing.
    precg_micro_bvftd <- mc_so@meta.data %>%
        as_tibble(rownames = "id")

    test_wilcox_tb <- wilcox_test_cols(precg_micro_bvftd, "clinical_dx", 
        list("age", "rin", "percent_mito", "pmi_h", "number_umi", "number_genes", "finalsite", "sex", "prep"))
 
    # boxplots of metadata, grouped by dx.
    meta_barplots <- ggplot_grouped_boxplot(precg_micro_bvftd, "clinical_dx",
        c("age", "rin", "percent_mito", "pmi_h", "number_umi", "number_genes", "prep"))
   
    # dump metadata 
    write_csv(precg_micro_bvftd, paste0(out_table, "meta.csv"))

    # dump dx
    write_csv(lm_tb_dx, paste0(out_table, "dx.csv"))

    # dump dx + umi
    write_csv(lm_tb_dx_umi, paste0(out_table, "dx_umi.csv"))

    # dump dx + mito + umi
    write_csv(lm_tb_dx_umi_mito, paste0(out_table, "dx_umi_mito.csv"))

    # dump full model
    write_csv(lm_tb, paste0(out_table, "lm_full.csv"))
    write_csv(test_wilcox_tb, paste0(out_table, "mann-whitney_metadata.csv"))
    
    # dump per-variable models
    walk(flatten(individual_vars), function(tb) {
        dx <- unique(tb$dx)
        term <- unique(tb$term)
        out_path <- str_glue("{out_table}_individual_{dx}_{term}.csv")
        writeLines(out_path)
        write_csv(tb, out_path)
    })

    # dump partial models
    write_csv(lm_pmi_age_site_sex_tb, str_glue("{out_table}_pmodel1.csv"))
    write_csv(lm_pmi_pct_mito_tb, str_glue("{out_table}_pmodel2.csv"))

    # rin + pmi models
    write_csv(lm_ctl_rin_umi_tb, str_glue("{out_table}_ctl_rin_umi.csv"))
    
    # write metadata barplots
    pdf(paste0(out_graph, "meta_barplots.pdf"), width = 16, height = 8 * (length(meta_barplots) / 2))
    tryCatch({plot(wrap_plots(meta_barplots, ncol = 2))}, error = print)
    dev.off()
}


run_lm_de <- function(
    seurat_obj,
    model_design,
    down_sample_cells = NULL,
    cores = NULL,
    ret_lm = NULL){
 
    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        seurat_obj <- subset(seurat_obj, cells = sample(rownames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clincal_dx, cluster_cell_type) %>% table %>% print
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
    
    if (!is.null(ret_lm) && ret_lm) {
        run_lm(expr_m, test_vars, model)
    } else {
        run_lm_broom(expr_m, test_vars, model)
    }
}

run_lm_broom <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr) {
        dat <- data.frame(expression = expr, test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lm(model, data = dat))
        },
        error = function(x) { NA })
    })
    return(lm_out_obj_l)
}

run_lm <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr) {
        dat <- data.frame(expression = expr, test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            lm(model, data = dat)
        },
        error = function(x) { NA })
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

nestedrank_test_cols <- function(tb, treatment_name, group_name, test_list) {
    lapply(test_list, function(var_name) {
        treatment_vec <- tb[[treatment_name]] %>%
            as.factor() %>% droplevels
        test_vec <- tb[[var_name]]
        if (!is.numeric(test_vec)) {
            test_vec <- as.numeric(factor(test_vec))
            var_name <- paste0("as_numeric_", var_name)
        }
        group_vec <- tb[[group_name]]

        out_tb <- broom::tidy(
            nestedRanksTest(x = treatment_vec, y = test_vec, group = group_vec))
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
