set.seed(27)
require(Seurat)
require(tidyverse)
require(future.apply)
require(lme4)
require(lmerTest)
require(broom.mixed)
options(future.globals.maxSize = Inf)

# in_seurat_expr <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_expression_matrix.rdat"
# in_seurat_meta <- "../analysis/seurat/20200613/tables/seurat_FTD_control_preCG_microglia_metadata.csv"
in_seurat_rds <-
  "../analysis/pci_import/pci_seurat.rds"

script_name <- "clinical_dx_de_microglia.R"
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
out_path_base <- "../analysis/clinical_dx_de_lme_microglia/"
out_table <- paste0(
  out_path_base, date,
  "/tables/clinical_dx_de_lme_microglia_")
dir.create(dirname(out_table), recursive = TRUE)
print(out_table)

main <- function() {
    nd_so <- readRDS(in_seurat_rds)

    cell_type <- "microglia"
    region <- "preCG"

    # Filter metadata to get cell ids for matching region, dx, cluster_cell_type.
    # Subset seurat object.
    cell_id_subset <- nd_so@meta.data %>%
        filter(cell_type == {{cell_type}},
               region == {{region}}) %>%
        pluck("cell_ids")
    mc_so <- subset(nd_so, cells = cell_id_subset)
    meta <- mc_so@meta.data
    rm(nd_so)
    gc()
    mc_so[["log_number_umi"]] <- log(mc_so[["number_umi"]])

    # Run model.
    plan(multicore, workers = 12)
    model_design <- "expression ~ clinical_dx + pmi + age + rin + sex + seq_batch + number_umi + percent_mito + Reads.Mapped.Antisense.to.Gene + Fraction.Reads.in.Cells + (1 | library_id)"
    system.time(lm_broom <- run_lmer_de(mc_so, model_design, down_sample_cells = 10000))
    print(sum(!is.na(lm_broom)))

    # Format output + filter by percent genes detected.
    plan(multicore)
    system.time({
        prop_detected <- prop_detected_generate(
            GetAssayData(mc_so, slot = "data"),
            mc_so@meta.data
        )
    })

    lm_broom_nna <- lm_broom[!is.na(lm_broom)]
    system.time({
        lm_tb_0 <- filter_and_format_lm_output(lm_broom_nna, mc_so, model_design,
            beta_regex = "clinical_dx", prop_detected, prop_detected_filter = 0)
    })
    system.time({
        lm_tb_5 <- filter_and_format_lm_output(lm_broom_nna, mc_so, model_design,
            beta_regex = "clinical_dx", prop_detected, prop_detected_filter = 0.05)})
    system.time({
        lm_tb_10 <- filter_and_format_lm_output(lm_broom_nna, mc_so, model_design,
            beta_regex = "clinical_dx", prop_detected, prop_detected_filter = 0.10)
    })
    plan(sequential)


    # Write out lm_broom. Format and write out lm_tb.
    saveRDS(lm_broom, paste0(out_table, "lm.rds"))
    saveRDS(lm_tb_0, paste0(out_table, "lm_0.rds"))
    saveRDS(lm_tb_5, paste0(out_table, "lm_5.rds"))
    saveRDS(lm_tb_10, paste0(out_table, "lm_10.rds"))

    write_csv(lm_tb_0, paste0(out_table, "lm_wider_tb_0.csv"))
    write_csv(lm_tb_5, paste0(out_table, "lm_wider_tb_5.csv"))
    write_csv(lm_tb_10, paste0(out_table, "lm_wider_tb_10.csv"))
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
    # If dx is present relevel seurat_obj so that control is the 1st factor level.
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


prop_cells_detected_dx <- function(expr, meta, clinical_dx_val) {
    is_dx <- meta[["clinical_dx"]] == clinical_dx_val
    prop_detected_dx <- future_apply(expr, 1, function(row) {
        row_dx <- row[is_dx]
        prop_detected <- (sum(row_dx > 0) / length(row_dx))
        prop_detected[!is.finite(prop_detected)] <- 0
        return(prop_detected)
    }) %>%
    round(3) %>%
    enframe(name = "gene", value = paste0("prop_detected_", clinical_dx_val))
    return(prop_detected_dx)
}

prop_detected_generate <- function(expr, meta) {
    # proportion of contrl cells expressing gene
    ctl_detect <- prop_cells_detected_dx(expr, meta, "Control")

    # proportion of dx cells expressing gene
    ad_detect <- prop_cells_detected_dx(expr, meta, "AD")
    ftd_detect <- prop_cells_detected_dx(expr, meta, "bvFTD")
    psp_detect <- prop_cells_detected_dx(expr, meta, "PSP-S")

    return(reduce(list(ad_detect, ftd_detect, psp_detect),
           inner_join,
           by = "gene",
           .init = ctl_detect))
}

filter_and_format_lm_output <- function(
    lm_out_obj_l,
    seurat_obj,
    model_design,
    beta_regex,
    prop_detected,
    prop_detected_filter = 0){

    # format lm output into tibble
    lm_tb <- lm_out_obj_l %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))

    # pivot dx output.
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = c("estimate", "p.value", "std.error", "statistic", "df"),
        names_glue = "{term}.{.value}"
    ) %>% arrange(term)
    lm_wd_tb <- pivot_wider_spec(lm_tb, id_cols = "gene", lm_wider_spec)

    if (nrow(lm_wd_tb) == 0)
        stop("lm output is of length zero. Check lm_out_obj_l and beta_regex arguments")

    meta <- seurat_obj@meta.data


    # join prop_detected with lm_wd_tb
    lm_out_prop <- inner_join(lm_wd_tb, prop_detected, by = "gene")

    # filter to genes expressed in greater than X prop of cells
    print(paste0("number of genes before prop detected filter: ", nrow(lm_out_prop)))
    lm_filter_tb <- lm_out_prop %>%
        filter(prop_detected_AD > prop_detected_filter |
               prop_detected_Control > prop_detected_filter |
               prop_detected_bvFTD > prop_detected_filter |
               `prop_detected_PSP-S` > prop_detected_filter )
    print(paste0("number of genes after prop detected filter: ", nrow(lm_filter_tb)))

    # fdr_pvalue correct
    lm_filter_fdr <- lm_filter_tb %>%
        mutate(
            clinical_dxAD.fdr.p.value = p.adjust(clinical_dxAD.p.value, method = "BH"),
            clinical_dxbvFTD.fdr.p.value = p.adjust(clinical_dxbvFTD.p.value, method = "BH"),
            `clinical_dxPSP-S.fdr.p.value` = p.adjust(`clinical_dxPSP-S.p.value`, method = "BH")
        )

    # print pval counts after correction
    lm_filter_fdr %>%
        summarize(across(contains("p.value"), ~sum(.x < 0.05))) %>%
        print(width = Inf)

    # formating
    lm_filter_out <- lm_filter_fdr %>%
        dplyr::select(
            gene,
            contains("estimate"),
            contains("p.value"),
            contains("statistic"),
            contains("prop_detected"),
            contains("fdr.p.value")
        ) %>%
        mutate(
            model = model_design,
            region = paste0(unique(meta$region)),
            cell_type = paste0(unique(meta$cell_type))
        )

    return(lm_filter_out)
}

if (interactive()) {
    main()
}
