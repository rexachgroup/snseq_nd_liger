# Compare one dx vs. all other conditions.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_subcluster_celltype_filter <- 
    "../../analysis/seurat_lchen/liger_subcluster_celltype_filter/celltype_filter.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_dx_isol/"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base, recursive = TRUE)

RESOURCES <- list(
    ncpus = 1,
    memory = 80,
    walltime = 172800,
    measure.memory = TRUE
)
prop_detected_filter <- 0.1
chunk_size <- 20

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    liger_meta <- readRDS(in_seurat_liger)
    subcluster_celltype_filter <- readRDS(in_subcluster_celltype_filter)

    # Filter metadata to get cell ids for matching region, cell_type.
    # Subset seurat object.
    liger_meta_subset <- liger_meta %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_subcluster = factor(ct_subcluster, levels = str_sort(unique(ct_subcluster), numeric = TRUE))
        ) %>%
        filter(ct_subcluster %in% subcluster_celltype_filter) %>%
        filter(cluster_cell_type  == "microglia", region == "preCG", liger_clusters == 7) %>%
        mutate(clinical_dx = fct_relevel(fct_other(clinical_dx, keep = "AD", other_level = "not_AD"), "not_AD"))
 

    # Split by subcluster + add model design. 
    model_designs <- c("expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    subcluster_wk <- liger_meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(region, cluster_cell_type, ct_subcluster, model_design) %>%
        group_nest(keep = TRUE)

    clearRegistry()
    batchExport(list(run_lmer_de = run_lmer_de, run_lmer = run_lmer,
                     prop_detected_generate = prop_detected_generate, prop_cells_detected_dx = prop_cells_detected_dx), reg = reg)
    ids <- batchMap(subcluster_worker,
        args = list(meta = subcluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    subcluster_wk$job_id <- getJobTable()$job.id
    subcluster_wk$prop_detected_filter <- prop_detected_filter
    submitJobs(ids, resources = RESOURCES)
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_wk.rds"))
    waitForJobs()

    # Fetch results.
    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join = pmap(list(data, model_design, job_id), function(data, model_design, job_id) {
                gc()
                broom_list <- loadResult(job_id)
                writeLines(paste(job_id))
                format_lm_output(broom_list, data, model_design, "clinical_dx")
            }))
    
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_wk.rds"))

    # Merge all subcluster results into one large table.
    subcluster_tb <- bind_rows(subcluster_wk$broom_join[map_lgl(subcluster_wk$broom_join, ~!any(is.na(.))) ])

    # Write out lm_broom lists.
    saveRDS(subcluster_tb, file.path(out_path_base, "subcluster_lme.rds"))
    write_csv(subcluster_tb, file.path(out_path_base, "subcluster_lme.csv"))
}

subcluster_worker <- function(meta, prop_detected_filter) {
    options(future.globals.maxSize = Inf)
    envdata <- new.env()
    load(file.path("..", unique(meta$filepath)), envdata)

    # Calculate proportion detected.
    prop_detected_tbl <- prop_detected_generate(meta, envdata$liger_so)
    detected <- prop_detected_tbl %>% filter(prop_detected_AD > prop_detected_filter |
                      prop_detected_Control > prop_detected_filter |
                      prop_detected_bvFTD > prop_detected_filter |
                      `prop_detected_PSP-S` > prop_detected_filter)

    # Subset.
    subcluster_so <- subset(envdata$liger_so, cells = meta$cell_ids, features = detected$gene)
    rm(envdata)
    gc()
    return(run_lmer_de(subcluster_so, meta, down_sample_cells = 10000, cores = RESOURCES$ncores))
}

run_lmer_de <- function(
    seurat_obj,
    meta,
    down_sample_cells = NULL,
    cores = NULL){

    model_design <- unique(meta$model_design)

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
    test_vars <- meta

    run_lmer(expr_m, test_vars, model, cores)
}

run_lmer <- function(expr_m, test_vars, model, cores = NULL) {
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr_r) {
        dat <- data.frame(expression = as.vector(expr_r), test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lmerTest::lmer(model, data = dat))
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
    writeLines(str_glue("prop_detected: {unique(meta$ct_subcluster)}, {paste0(dim(subcluster_so), collapse = ' ')}"))

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
    liger_meta,
    model_design,
    beta_regex){

    # format lm output into tibble
    lm_tb <- lm_out_obj_l[!is.na(lm_out_obj_l)] %>%
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

if (!interactive()) {
    main()
}
