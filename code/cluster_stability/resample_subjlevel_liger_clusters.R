# Resample + recluster at subject id level. Sample 80% of samples per dx.

liblist <- c("tidyverse", "Seurat", "liger", "batchtools")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

# in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
subcluster_cts <- "../../analysis/liger_subcluster/20200824/"

out_path_base <- "../../analysis/seurat_lchen/cluster_stability_subjlevel/"
dir.create(out_path_base)
batchtools <- file.path(out_path_base, "batchtools")

RESOURCES <- list(
    ncpus = 16,
    memory = 64,
    walltime = 36000,
    measure.memory = TRUE
)

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, writeable = TRUE)
        removeRegistry()
    } 
    reg <- makeRegistry(batchtools, seed = 1)

    meta <- as_tibble(readRDS(in_seurat_meta))

    meta_subset <- meta %>% 
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-")) %>%
        filter((region == "preCG" & cluster_cell_type == "microglia") | (region == "insula" & cluster_cell_type == "excitatory"))
    
    subcluster_tbs <- meta_subset %>%
        group_by(region, cluster_cell_type) %>%
        group_nest(keep = TRUE)

    subcluster_wk <- subcluster_tbs %>%
        mutate(resample_cellids = map(data, function(data) {
            resample_subjdx(data, 5)
        })) %>%
        unnest(resample_cellids)
    
    clearRegistry()
    ids <- batchMap(liger_cluster_sobj, 
        args = list(
            cell_ids = subcluster_wk$resample_cellids,
            data = subcluster_wk$data
        ), 
        more.args = list(
            liblist = liblist
        ))
    subcluster_wk$job.id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    waitForJobs()
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_resamples.rds"))
}
 
resample_subjdx <- function(meta, ntests = 5, nprop = 0.8, cell_id_col = "cell_ids") {
    set.seed(1)
    return(map(1:ntests, function(i) {
        autopsy_ids <- meta %>%
            group_by(autopsy_id) %>%
            summarize(autopsy_id = unique(autopsy_id), clinical_dx = unique(clinical_dx))
        sampled_ids <- autopsy_ids %>%
            group_by(clinical_dx) %>%
            slice_sample(prop = nprop)
        cell_ids <- meta %>%
            filter(autopsy_id %in% sampled_ids$autopsy_id) %>%
            pluck(cell_id_col)
        return(cell_ids)
    }))
}

# Worker function: subset sobj to this instance's cell_ids, then convert to liger object and run clustering
liger_cluster_sobj <- function(cell_ids, data, liblist) {
    filepath <- unique(data$filepath)
    k <- 20
    lambda <- 1

    stopifnot("Duplicate cell ids cause errors in seuratToLiger -- don't use bootstrapping"=any(!duplicated(cell_ids)))
    l <- lapply(liblist, require, character.only = TRUE)

    envdata <- new.env()
    load(file.path("..", filepath), envdata) 
    subset_sobj <- subset(envdata$liger_so, cells = cell_ids)

    lo <- seuratToLiger(subset_sobj,
        combined.seurat = TRUE, names = "use-meta",
        meta.var = "prep", assays.use = NULL, raw.assay = "RNA",
        remove.missing = FALSE, renormalize = FALSE, use.seurat.genes = FALSE,
        num.hvg.info = NULL, use.idents = FALSE, use.tsne = TRUE, cca.to.H = FALSE) 

    lo <- lo %>%
        normalize %>%
        selectGenes %>%
        scaleNotCenter %>%
        optimizeALS(k = k, lambda = lambda, thresh = 5e-5, nrep = 3) %>%
        quantileAlignSNF %>%
        runUMAP

    return(lo)
}

if (!interactive()) {
    main()
}
