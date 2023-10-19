liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

META <- "../../analysis/seurat_lchen/"
SUBCLUSTER_FILTERED <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/01_subcluster_filtered_meta.csv"
SOBJ <- "../../analysis/pci_import/pci_seurat.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_excluded_dge/subcluster_dge/"

batchtools <- file.path(OUT_DIR, "batchtools")

RESOURCES <- list(
    ncpus = 1,
    memory = 80,
    walltime = 172800,
    measure.memory = TRUE,
    chunks.as.arrayjobs = TRUE,
    partition = "bigmem"
)

prop_detected_filter <- 0.1
model_designs <- c("expression ~ clinical_dx + pmi + age + sex + number_umi + percent_mito + (1 | library_id)")
chunk_size <- 5
main <- function() {
    dir.create(OUT_DIR, recursive = TRUE)
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, conf.file = "../batchtools.conf.R", writeable = TRUE)
        #removeRegistry(reg)
        #sweepRegistry()
    } else {
        reg <- makeRegistry(batchtools, packages = liblist, conf.file = "../batchtools.conf.R")
    }
    source("../lib/subcluster_dge.R")
    liger_meta <- read_csv(SUBCLUSTER_FILTERED)

    liger_meta_subset <- liger_meta %>%
        filter(
            cluster_cell_type == "astrocyte",
            clinical_dx %in% c("PSP-S", "Control"),
        ) %>%
        mutate(
            clinical_dx = fct_relevel(clinical_dx, "Control"),
            library_id = as.factor(library_id)
        )
    
    subcluster_wk <- liger_meta_subset %>%
        cross_join(tibble(model_design = model_designs)) %>%
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
    saveRDS(subcluster_wk, file.path(OUT_DIR, "subcluster_wk.rds"))
    waitForJobs()

    # Fetch results.
    subcluster_tb <- subcluster_wk %>%
        filter(job_id %in% findDone()$job.id) %>%
        mutate(broom_join = pmap(list(data, model_design, job_id), function(data, model_design, job_id) {
                gc()
                broom_list <- loadResult(job_id)
                writeLines(paste(job_id))
                if (all(is.na(broom_list))) {
                    return(NA)
                }
                format_lm_output(broom_list, data, model_design, "clinical_dx")
            }))
    
    saveRDS(subcluster_tb, file.path(OUT_DIR, "subcluster_wk.rds"), compress = FALSE)

    # Merge all subcluster results into one large table and save.
    subcluster_tbs <- bind_rows(subcluster_tb$broom_join[map_lgl(subcluster_tb$broom_join, ~!any(is.na(.))) ])

    saveRDS(subcluster_tbs, file.path(OUT_DIR, "subcluster_lme.rds"), compress = FALSE)

    # For csv export, write each subcluster to a separate file.
    subcluster_tb <- subcluster_tb %>%
        mutate(filepath = file.path(
            OUT_DIR, 
            "lme_tables",
            paste(region, cluster_cell_type, sep = "-"),
            paste0(ct_subcluster, ".csv"))
        )

    pwalk(list(subcluster_tb$filepath, subcluster_tb$broom_join), function(filepath, tb) {
            dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
            if (!all(is.na(tb))) {
                write_csv(tb, filepath)
            }
        })
    
    subcluster_gene_summary <- subcluster_tb %>%
        filter(!is.na(broom_join)) %>%
        mutate(gene_cts = pmap(., function(...) {
                cr <- list(...)
                summarize_gene_counts(cr$broom_join)
            })) %>%
        select(region, cluster_cell_type, ct_subcluster, gene_cts)

    subcluster_gene_summary %>%
        unnest(gene_cts, names_repair = "universal") %>%
        select(region, cluster_cell_type, ct_subcluster = ct_subcluster...3, dx, up, down) %>%
        write_csv(file.path(OUT_DIR, "dx_dge_summary.csv"))

    
}

if (!interactive()) main()
