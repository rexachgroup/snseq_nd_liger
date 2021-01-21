# clustering on liger clusters using VariancePartition / dream.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_dream/"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base)

RESOURCES <- list(
    ncpus = 16,
    memory = 64,
    walltime = 36000,
    measure.memory = TRUE
)

main <- function() {
    if (dir.exists(batchtools)) {
        loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        removeRegistry(wait = 0)
    }
    reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    liger_meta <- readRDS(in_seurat_liger)

    # Filter metadata to get cell ids for matching region, cell_type.
    # Subset seurat object.
    liger_meta_subset <- liger_meta %>%
        mutate(
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        filter(
            (cluster_cell_type == "microglia" & region == "preCG"),
            liger_clusters %in% c(1, 7)
        )

    # Split by subcluster + add model design. 
    model_designs <- c("~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")

    subcluster_tbs <- liger_meta_subset %>%
        inner_join(tibble(model_design = model_designs), by = character()) %>%
        arrange(ct_subcluster, model_design) %>%
        group_by(ct_subcluster, model_design) %>%
        group_nest(keep = TRUE)

    # Filter by proportion detected per dx.
    prop_detected_filter <- 0.1
    subcluster_wk <- subcluster_tbs %>%
        mutate(prop_detected = map(data, function(data) {
            envdata <- new.env()
            load(file.path("..", unique(data$filepath)), envdata) 
            subcluster_so <- subset(envdata$liger_so, cells = data$cell_ids)
            writeLines(str_glue("prop_detected: {unique(data$ct_subcluster)}"))

            prop_detected <- prop_detected_generate(
                GetAssayData(subcluster_so, slot = "data"), data
            )

            prop_detected %>%
                filter(prop_detected_AD > prop_detected_filter |
                       prop_detected_Control > prop_detected_filter |
                       prop_detected_bvFTD > prop_detected_filter |
                       `prop_detected_PSP-S` > prop_detected_filter )
        })) 
    
    # Submit DREAM runs to cluster nodes.
    clearRegistry()
    batchExport(list(run_dream_de = run_dream_de))
    ids <- batchMap(subcluster_worker, subcluster_wk$data, subcluster_wk$prop_detected)
    subcluster_wk$job_id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_dream_wk.rds"))
    waitForJobs()

    # Extract output from cluster results.
    # Use topTable() to get results for each contrast.
    subcluster_wk <- subcluster_wk %>%
        mutate(broom_join = pmap(list(data, job_id), function(data, job_id) {
                print(unique(data$ct_subcluster))
                dream_obj <- loadResult(job_id)
                dx_coefs <- c("clinical_dxbvFTD", "clinical_dxAD", "clinical_dxPSP-S")

                dx_coefs <- setNames(dx_coefs, dx_coefs)
                dx_tbls <- map(dx_coefs, function(x) {
                    topTable(dream_obj, x, number = nrow(dream_obj)) %>%
                        as_tibble(rownames = "gene")
                }) %>%
                bind_rows(.id = "dx") 
            })
        )

    # Merge all subcluster results into one large table.
    subcluster_tb <- subcluster_wk %>% 
        select(ct_subcluster, model_design, broom_join) %>%
        unnest(broom_join) %>%
        pivot_wider(names_from = "dx", values_from = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "z.std"))

    # Write out lists.
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_dream_objs.rds"))
    write_csv(subcluster_tb, file.path(out_path_base, "subcluster_dream.csv"))
}

subcluster_worker <- function(meta, prop_detected) {
    envmeta <- new.env()
    load(file.path("..", unique(meta$filepath)), envmeta) 
    subcluster_so <- subset(envmeta$liger_so, cells = meta$cell_ids, features = prop_detected$gene)
    return(run_dream_de(subcluster_so, unique(meta$model_design), down_sample_cells = 10000))
}


run_dream_de <- function(seurat_obj, model_design, down_sample_cells = NULL, ncores = 16){
    # Downsample if down_sample_cells defined.
    if (!is.null(down_sample_cells) && ncol(seurat_obj) > down_sample_cells) {
        writeLines(str_glue("downsampling to {down_sample_cells}"))
       seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), down_sample_cells))
        seurat_obj@meta.data %>%
            select(region, clinical_dx, cluster_cell_type) %>% table %>% print
    }

    expr_m <- GetAssayData(seurat_obj, slot = "counts")

    test_vars <- as.data.frame(seurat_obj@meta.data)

    # Create contrasts for clinical_dx.
    rownames(test_vars) <- seurat_obj@meta.data$cell_ids

    dgelist <- DGEList(expr_m)
    dgelist_norm <- calcNormFactors(dgelist)
    
    # browser()
    if ("clinical_dx" %in% colnames(test_vars)) {
        test_vars$clinical_dx <- fct_relevel(test_vars$clinical_dx, "Control")
        #         contrastAD <- getContrast(expr_m, model_design, test_vars, c("clinical_dxAD"))
        #         contrastbvFTD <- getContrast(expr_m, model_design, test_vars, c("clinical_dxbvFTD"))
        #         contrastPSP.S <- getContrast(expr_m, model_design, test_vars, c("clinical_dxPSP-S"))
        #         contrast_list <- cbind(contrastAD, contrastbvFTD, contrastPSP.S)
    }

    elist_voom <- voomWithDreamWeights(counts = dgelist_norm, formula = model_design, data = test_vars, BPPARAM = SnowParam(ncores))
    dx_fit <- dream(elist_voom, model_design, test_vars, BPPARAM = SnowParam(ncores), ddf = "Kenward-Roger")

    return(dx_fit)
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

if (!interactive()) {
    main()
}
