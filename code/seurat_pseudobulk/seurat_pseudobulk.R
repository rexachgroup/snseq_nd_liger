set.seed(0)
liblist <- c("Seurat", "SeuratObject", "tidyverse", "batchtools", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_subcluster_celltype_filter <- 
    "../../analysis/seurat_lchen/liger_subcluster_celltype_filter/celltype_filter.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/seurat_pseudobulk/"
batchtools <- file.path(out_path_base, "batchtools")
pbulk_out <- file.path(out_path_base, "pbulk_sobj.rds")
dir.create(out_path_base, recursive = TRUE)

#prop_detected_filter <- 0.1
#chunk_size <- 40
    
model_design <- "expression ~ clinical_dx + pmi + age + sex + mean_number_umi + mean_percent_mito"

RESOURCES <- list(
    ncpus = 1,
    memory = 32,
    walltime = 172800,
    measure.memory = TRUE
)

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    reg$packages <- liblist
    sobj <- readRDS(in_seurat_rds)
    pbulk_sobj <- AggregateExpression(sobj, return.seurat = T, group.by = c("library_id", "cluster_cell_type"))
    saveRDS(pbulk_sobj, pbulk_out, compress = F)

    meta <- sobj@meta.data %>%
        group_by(library_id, cluster_cell_type) %>%
        summarize(
            clinical_dx = unique(clinical_dx),
            pmi = unique(pmi),
            age = unique(age),
            sex = unique(sex),
            mean_number_umi = mean(number_umi),
            mean_percent_mito = mean(percent_mito),
            .groups = "drop"
        )
    meta_bind <- meta %>% mutate(ct_lib = str_glue("{library_id}_{cluster_cell_type}"))
    pbulk_meta <- pbulk_sobj@meta.data %>%
        rownames_to_column("ct_lib") %>%
        left_join(meta_bind, by = "ct_lib") %>%
        column_to_rownames("ct_lib")
    pbulk_sobj@meta.data <- pbulk_meta
    saveRDS(pbulk_sobj, pbulk_out, compress = F)

    cluster_wk <- pbulk_meta %>%
        mutate(
            filepath = pbulk_out,
            model_design = model_design,
            cell_ids = rownames(pbulk_meta)
        ) %>%
        group_by(cluster_cell_type, model_design) %>%
        group_nest(keep = T)
    
    clearRegistry()
    batchExport(list(run_lm_de = run_lm_de, run_lm_broom = run_lm_broom), reg = reg)
    ids <- batchMap(cluster_worker,
        args = list(meta = cluster_wk$data)
    )
    ids <- findNotDone()
    cluster_wk$job_id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    saveRDS(cluster_wk, file.path(out_path_base, "cluster_wk.rds"))
    waitForJobs()

    cluster_wk <- cluster_wk %>%
        mutate(dge_long = map(job_id, function(i) {
            res <- loadResult(i)
            fmt_dge_long(res, "clinical_dx")
        }))

    cluster_wk <- cluster_wk %>%
        mutate(dge_wide = map(dge_long, fmt_dge_wide, "clinical_dx"))
    
    saveRDS(cluster_wk, file.path(out_path_base, "cluster_wk.rds", compress = F))

    pwalk(cluster_wk, function(...) {
        cr <- list(...)
        out_l_path <- str_glue("{out_path_base}/{cr$cluster_cell_type}_long.csv")
        out_w_path <- str_glue("{out_path_base}/{cr$cluster_cell_type}_wide.csv")
        write_csv(cr$dge_long, out_l_path)
        write_csv(cr$dge_wide, out_w_path)
    })
    
}

cluster_worker <- function(meta) {
    options(future.globals.maxSize = Inf)
    sobj <- readRDS(file.path(unique(meta$filepath)))
    
    # Subset.
    subset_sobj <- subset(sobj, cells = meta$cell_ids)

    return(run_lm_de(subset_sobj, unique(meta$model_design), down_sample_cells = NULL, cores = RESOURCES$ncores))
}


run_lm_de <- function(
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

    run_lm_broom(expr_m, test_vars, model, cores)
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

fmt_dge_long <- function(dge_list, dx_regex) {
    bind_tb <- bind_rows(dge_list, .id = "gene_name")
    dx_cols <- str_subset(unique(bind_tb$term), dx_regex)
    p_adjust_dx <- lapply(dx_cols, function(col) {
        bind_tb %>%
            filter(term == col) %>%
            mutate(p.val.adj = p.adjust(p.value, method = "fdr")) %>%
            select(gene_name, term, p.val.adj)
    }) %>% bind_rows()
    bind_tb <- left_join(bind_tb, p_adjust_dx, by = c("gene_name", "term"))
    return(bind_tb)
}

fmt_dge_wide <- function(dge_tb_long, dx_regex) {
    dge_tb_long %>%
        filter(str_detect(term, dx_regex)) %>%
        pivot_wider(names_from = "term", values_from = c("estimate", "std.error", "statistic", "p.value", "p.val.adj"))
}

main()
