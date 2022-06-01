# Export aveExpr / proportion genes per library for liger microglia clusters.
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_liger_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"

out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_counts_export"
batchtools <- file.path(out_path_base, "batchtools")
dir.create(out_path_base, recursive = TRUE)

RESOURCES <- list(
    ncpus = 1,
    memory = 32,
    walltime = 3600,
    measure.memory = TRUE
)
chunk_size <- 30

main <- function(){
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
        reg$packages <- liblist
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "batchtools.conf.R")
    }
    liger_meta <- readRDS(in_liger_meta)

    subcluster_wk <- liger_meta %>%
        filter(
            region == "preCG",
            cluster_cell_type == "microglia",
            liger_clusters %in% c("0", "1", "4", "7")
        ) %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(ct_subcluster, library_id) %>% group_nest(keep = TRUE)

    clearRegistry()
    batchExport(mget(ls()))

    ids <- batchMap(subcluster_worker, args = list(meta = subcluster_wk$data))
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    subcluster_wk$job.id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    waitForJobs()

    subcluster_wk <- subcluster_wk %>%
        mutate(gene_tb = map(job.id, function(x) {
            res <- loadResult(x)
            inner_join(res[[1]], res[[2]], by = "gene") %>%
                rename(prop_detected_lib = prop_detected_all, aveExpr_lib = means)
        }))

    subcluster_wk %>%
        select(-data, -job.id) %>%
        unnest(gene_tb) %>%
        group_split(ct_subcluster) %>%
        walk(function(x) {
            out_path <- file.path(out_path_base, paste0(unique(x$ct_subcluster), "-library_counts.csv"))
            writeLines(out_path)
            print(x)
            write_csv(x, out_path)
        })
}

subcluster_worker <- function(meta) {
    # Load in subcluster file.
    options(future.globals.maxSize = Inf)
    envdata <- new.env()
    load(file.path(unique(meta$filepath)), envdata)
    
    # Subset.
    subcluster_so <- subset(envdata$liger_so, cells = meta$cell_ids)
    rm(envdata)
    gc()

    print(dim(subcluster_so))
    prop_detected_tbl <- prop_detected_generate(meta, subcluster_so)
    ave_expr_tbl <- ave_expr_generate(meta, subcluster_so)
    return(list(prop_detected_tbl, ave_expr_tbl))
}

prop_detected_generate <- function(meta, subcluster_so) {
    writeLines(str_glue("prop_detected: {paste0(unique(meta$ct_subcluster), collapse = ' ')}, {paste0(dim(subcluster_so), collapse = ' ')}"))
    expr_m <- GetAssayData(subcluster_so, slot = "data")

    expr_cells_gz <- rowSums(expr_m > 0)
    frac_gz <- as.double(expr_cells_gz) / ncol(expr_m)
    prop_detected_all <- tibble(gene = rownames(expr_m), prop_detected_all = frac_gz)

    return(prop_detected_all)
}

ave_expr_generate <- function(meta, subcluster_so) {
    writeLines(str_glue("ave_expr_generate: {paste0(unique(meta$ct_subcluster), collapse = ' ')}, {paste0(dim(subcluster_so), collapse = ' ')}"))
    expr_m <- GetAssayData(subcluster_so, slot = "data")
    means <- rowMeans(expr_m)
    return(tibble(gene = rownames(expr_m), means = means))
}

if (!interactive()) {
    main()
}
