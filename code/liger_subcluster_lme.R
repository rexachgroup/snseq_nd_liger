ssset.seed(27)
require(Seurat)
require(tidyverse)
require(future.apply)
require(lme4)
require(lmerTest)
require(broom.mixed)
options(future.globals.maxSize = Inf)

in_seurat_rds <-
  "../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../analysis/seurat_lchen/liger_subcluster_metadata.rds"

script_name <- "clinical_dx_de_microglia.R"
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
out_path_base <- "../analysis/seurat_lchen/liger_subcluster_lme/"
out_table <- file.path(out_path_base)
dir.create(dirname(out_table), recursive = TRUE)
print(out_table)

main <- function() {
    # nd_so <- readRDS(in_seurat_rds)
    liger_meta <- readRDS(in_seurat_liger)

    #cell_type <- "microglia"
    #region <- "preCG"

    # Filter metadata to get cell ids for matching region, cell_type.
    # Subset seurat object.
    liger_meta_subset <- liger_meta %>%
        mutate(
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        filter(
            (cluster_cell_type == "microglia" & region == "preCG") | 
            (cluster_cell_type == "microglia" & region == "insula") |
            (cluster_cell_type == "excitatory" & region == "insula")
        )

    # Split by subcluster, then run lm_broom. 
    model_designs <- c("expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    subcluster_tbs <- liger_meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(ct_subcluster, model_design) %>%
        group_nest(keep = TRUE)

    # Run lme_broom.
    subcluster_wk <- subcluster_tbs %>%
        mutate(broom_lists = map(data, function(data) {
            nd_so <- readRDS(in_seurat_rds)
            subcluster_so <- subset(nd_so, cells = data$cell_ids)
            rm(nd_so)
            gc()
            
            nworkers = max(floor(80 / (nrow(subcluster_so) * 0.015)), 1)
            plan(multicore, workers = nworkers)
            system.time(lme_broom <- run_lmer_de(subcluster_so, data$model_design, down_sample_cells = 10000))
            return(lme_broom)
        }))

    # Format output + filter by proportion of genes detected.
    plan(multicore, workers = 8)

    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join = pmap(list(data, broom_lists), function(data, lme_broom) {
                tryCatch({
                    if (any(!is.na(lme_broom))) {
                        nd_so <- readRDS(in_seurat_rds)
                        subcluster_so <- subset(nd_so, cells = data$cell_ids)
                        rm(nd_so)
                        gc()

                        system.time({
                            prop_detected <- prop_detected_generate(
                                GetAssayData(subcluster_so, slot = "data"),
                                data
                            )
                        })

                        lme_broom_nna <- lme_broom[!is.na(lme_broom)]
                        print(unique(data$ct_subcluster))
                        system.time({
                            lme_tb_5 <- filter_and_format_lm_output(lme_broom_nna, subcluster_so, data, unique(data$model_design),
                                beta_regex = "clinical_dx", prop_detected, prop_detected_filter = 0.05)})
                        return(lme_tb_5)
                    }
                }, error = function(x) {print(x); return(NA)})
            })
        )

    plan(sequential)

    walk(subcluster_wk$broom_join, function(lme_tb) {
        if (!any(is.na(lme_tb))) {
            print(unique(lme_tb$ct_subcluster))
        }
    })

    # Merge all subcluster results into one large table.
    subcluster_tb <- bind_rows(subcluster_wk$broom_join[map_lgl(subcluster_wk$broom_join, ~!any(is.na(.))) ])

    # Write out lm_broom lists.
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_lme_objs.rds"))
    write_csv(subcluster_tb, file.path(out_path_base, "subcluster_lme.csv"))

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
        error = function(x) { NA })
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
    liger_meta,
    model_design,
    beta_regex,
    prop_detected,
    prop_detected_filter = 0){

    # format lm output into tibble
    lm_tb <- lm_out_obj_l %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))
    
    if (ncol(lm_tb) != 9) {
        browser()
    }

    if (nrow(lm_tb) == 0) {
        stop("lm output is of length zero. Check lm_out_obj_l and beta_regex arguments")
    }

    # pivot dx output.
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = c("estimate", "p.value", "std.error", "statistic", "df"),
        names_glue = "{term}.{.value}"
    ) %>% arrange(term)
    lm_wd_tb <- pivot_wider_spec(lm_tb, id_cols = "gene", lm_wider_spec)


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
            region = paste0(unique(liger_meta$region)),
            cell_type = paste0(unique(liger_meta$cluster_cell_type)),
            ct_subcluster = paste0(unique(liger_meta$ct_subcluster))
        )

    return(lm_filter_out)
}

if (interactive()) {
    main()
}
