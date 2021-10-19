# Compare excitatory layers 2 3 4 vs. 5 6.
set.seed(0)
options(deparse.max.lines = 5)
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <-
  "../../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/seurat_excitatory_layer_lme/"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base, recursive = TRUE, showWarnings = FALSE)

RESOURCES <- list(
    ncpus = 1,
    memory = 80,
    walltime = 172800,
    measure.memory = TRUE
)
chunk_size <- 1
prop_detected_filter <- 0.1
#model_design <- "expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)"
model_design <- "expression ~ cluster_celltype_layer"

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    reg$packages <- liblist
    meta <- readRDS(in_seurat_meta)

    # Filter metadata to get cell ids for matching region, cell_type.
    # Subset seurat object.
    meta_subset <- meta %>%
        mutate( 
            ct_cluster = paste(region, cluster_celltype_layer, sep = "-"),
            log_number_umi = log(number_umi),
            filepath = in_seurat_rds,
            model_design = model_design
        ) %>%
        mutate(
            ct_cluster = factor(ct_cluster, levels = str_sort(unique(ct_cluster), numeric = TRUE))
        ) %>%
        filter(cluster_celltype_layer %in% c("G234", "G56"))

    # Split by cluster, then run lm_broom. 
    cluster_wk <- meta_subset %>%
        group_by(region, clinical_dx) %>%
        group_nest(keep = TRUE)

    clearRegistry()
    batchExport(list(run_lm_de = run_lm_de, run_lm = run_lm,
                     prop_detected_generate = prop_detected_generate, prop_cells_detected_dx = prop_cells_detected_dx), reg = reg)
    ids <- batchMap(cluster_worker,
        args = list(meta = cluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    cluster_wk$job_id <- getJobTable()$job.id
    cluster_wk$prop_detected_filter <- prop_detected_filter
    submitJobs(ids, resources = RESOURCES)
    saveRDS(cluster_wk, file.path(out_path_base, "cluster_wk.rds"))
    waitForJobs()

    cluster_wk <- cluster_wk %>%
        mutate(broom_join = pmap(., function(...) {
                gc()
                current_row <- list(...)
                tryCatch({
                    broom_list <- loadResult(current_row$job_id)
                    writeLines(paste(current_row$job_id))
                    format_lm_output(broom_list, current_row$data, "cluster_celltype_layer")
                }, error = print)
            }))

    cluster_wk <- cluster_wk %>%
        mutate(filepath = file.path(
                out_path_base,
                "lm_tables",
                str_glue("layer_contrast_{region}_{clinical_dx}.csv")
            ))
    
    saveRDS(cluster_wk, file.path(out_path_base, "cluster_wk.rds"))

    pwalk(cluster_wk, function(...) {
        current_row <- list(...)
        dir.create(dirname(current_row$filepath), recursive = TRUE, showWarnings = FALSE)
        writeLines(current_row$filepath)
        if (is.data.frame(current_row$broom_join)) {
            write_csv(current_row$broom_join, current_row$filepath)
        }
    })

}

cluster_worker <- function(meta, prop_detected_filter) {
    options(future.globals.maxSize = Inf)
    cluster_so <- readRDS(file.path(unique(meta$filepath))) 
    
    # Subset.
    cluster_so <- subset(cluster_so, cells = meta$cell_ids)
    gc()
    
    prop_detected_tbl <- prop_detected_generate(meta, cluster_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_AD > prop_detected_filter |
                      prop_detected_Control > prop_detected_filter |
                      prop_detected_bvFTD > prop_detected_filter |
                      `prop_detected_PSP-S` > prop_detected_filter)
    cluster_so <- subset(cluster_so, cells = meta$cell_ids, features = detected$gene)

    return(run_lm_de(cluster_so, meta, unique(meta$model_design), cores = RESOURCES$ncores))
}

run_lm_de <- function(
    seurat_obj,
    meta,
    model_design,
    down_sample_cells = NULL,
    cores = NULL){

    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) 
    {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
        seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "data")

    # Convert model to formula.
    # If dx is present relevel seurat_obj so that control is the 1st factor level.
    model <- as.formula(model_design)
    test_vars <- meta

    if ("clinical_dx" %in% colnames(test_vars)) {
        test_vars <- mutate(test_vars, clinical_dx = fct_relevel(clinical_dx, "Control"))
    }

    run_lm(expr_m, test_vars, model, cores)
}


run_lm <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr_r) {
        dat <- data.frame(expression = as.vector(expr_r), test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lm(model, data = dat))
        },
        error = function(x) { print(x); return(NA) })
    })
    return(lm_out_obj_l)
}


prop_cells_detected_dx <- function(expr, meta, clinical_dx_val) {
    is_dx <- meta[["clinical_dx"]] == clinical_dx_val
    expr_genes_dx <- expr[,is_dx]
    expr_isexpressed <- rowSums(expr > 0)
    prop_detected <- expr_isexpressed / sum(is_dx)
    prop_detected_tb <- enframe(prop_detected, name = "gene", value = paste0("prop_detected_", clinical_dx_val))
    return(prop_detected_tb)
}

prop_detected_generate <- function(meta, subcluster_so) {
    expr <- GetAssayData(subcluster_so, slot = "data")

    # proportion of control cells expressing gene
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

format_lm_output <- function(
    lm_out_obj_l,
    meta,
    beta_regex){

    # format lm output into tibble
    lm_l <- lm_out_obj_l[!is.na(lm_out_obj_l)] 
    if (length(lm_l) == 0) {
        return(NA)
    }
    lm_tb <- lm_l %>%
        bind_rows(.id = "gene") %>%
        filter(grepl(beta_regex, term) | is.na(term))
    
    if (nrow(lm_tb) == 0) {
        stop("lm output is of length zero. Check lm_out_obj_l and beta_regex arguments")
    }

    # pivot dx output.
    lm_wider_spec <- build_wider_spec(lm_tb,
        names_from = "term",
        values_from = c("estimate", "p.value", "std.error", "statistic"),
        names_glue = "{term}.{.value}"
    ) %>% arrange(term)
    lm_wd_tb <- pivot_wider_spec(lm_tb, id_cols = "gene", lm_wider_spec)

    # fdr_pvalue correct
    lm_filter_fdr <- lm_wd_tb %>%
        mutate(across(contains("p.value"), ~p.adjust(.x, method = "BH"), .names = "{.col}.adj"))

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
            contains("statistic")
        ) %>%
        mutate(
            model = paste0(unique(meta$model_design), collapse = ":"),
            region = paste0(unique(meta$region), collapse = ":"),
            cell_type = paste0(unique(meta$cluster_cell_type), collapse = ":"),
            ct_cluster = paste0(unique(meta$ct_cluster), collapse=":")
        )

    return(lm_filter_out)
}

if (!interactive()) {
    main()
}
