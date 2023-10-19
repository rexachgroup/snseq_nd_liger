liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

META <- "../../analysis/seurat_lchen/"
SUBCLUSTER_FILTERED <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/01_subcluster_filtered_meta.csv"
SOBJ <- "../../analysis/pci_import/pci_seurat.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_excluded_dge/cluster_dge/"
batchtools <- file.path(OUT_DIR, "batchtools")

RESOURCES <- list(
    ncpus = 1,
    memory = 80,
    walltime = 172800,
    measure.memory = TRUE
)

prop_detected_filter <- 0.1
chunk_size <- 5
model_design <- "expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)"

main <- function() {
    dir.create(OUT_DIR) 
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    source("../lib/cluster_dge.R")
    reg$packages <- liblist
    meta <- read_csv(SUBCLUSTER_FILTERED)
    
    meta_subset <- meta %>%
        filter(
            cluster_cell_type == "astrocyte",
            clinical_dx %in% c("PSP-S", "Control"),
        ) %>%
        mutate(
            clinical_dx = fct_relevel(clinical_dx, "Control"),
            ct_cluster = paste(region, cluster_cell_type, sep = "-"),
            library_id = as.factor(library_id),
            model_design = model_design,
            filepath = SOBJ
        )
    
    cluster_wk <- meta_subset %>%
        group_by(region, cluster_cell_type, ct_cluster) %>%
        group_nest(keep = TRUE)
    
    clearRegistry()
    batchExport(list(run_lmer_de = run_lmer_de, run_lmer = run_lmer,
                     prop_detected_generate = prop_detected_generate, prop_cells_detected_dx = prop_cells_detected_dx), reg = reg)
    ids <- batchMap(cluster_worker,
        args = list(meta = cluster_wk$data),
        more.args = list(prop_detected_filter = prop_detected_filter)
    )
    ids <- findNotDone() %>%
        mutate(chunk = chunk(job.id, chunk.size = chunk_size))
    jt_ids <- getJobTable() %>%
        mutate(ct_cluster = map(job.pars, function(x) { unique(x$meta$ct_cluster) }) %>% unlist)
    cluster_wk <- left_join(cluster_wk, jt_ids, by = c("ct_cluster"))
    cluster_wk$prop_detected_filter <- prop_detected_filter
    submitJobs(ids, resources = RESOURCES)
    saveRDS(cluster_wk, file.path(OUT_DIR, "cluster_wk.rds"))
    waitForJobs()
    
    cluster_wk <- cluster_wk %>%
        filter(!is.na(job.id)) %>%
        mutate(broom_join = pmap(., function(...) {
                gc()
                current_row <- list(...)
                broom_list <- loadResult(current_row$job.id)
                writeLines(paste(current_row$job.id))
                format_lm_output(broom_list, current_row$data, "clinical_dx")
            }))

    cluster_wk <- cluster_wk %>%
        mutate(filepath = file.path(
                OUT_DIR,
                "lme_tables",
                paste0(ct_cluster, ".csv")
            ))
    
    saveRDS(cluster_wk, file.path(OUT_DIR, "cluster_wk.rds"))

    pwalk(cluster_wk, function(...) {
        current_row <- list(...)
        dir.create(dirname(current_row$filepath), recursive = TRUE, showWarnings = FALSE)
        writeLines(current_row$filepath)
        if (is.data.frame(current_row$broom_join)) {
            write_csv(current_row$broom_join, current_row$filepath)
        }
    })
}


