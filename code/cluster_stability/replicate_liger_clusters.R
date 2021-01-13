# replication of original liger clustering.

liblist <- c("tidyverse", "Seurat", "liger", "batchtools")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
subcluster_cts_rdata <- "../../analysis/liger_subcluster/20200824/"
in_seurat_rdata <- "../../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"

out_path_base <- "../../analysis/seurat_lchen/cluster_repl"
batchtools <- "../../analysis/seurat_lchen/cluster_repl/batchtools"

RESOURCES <- list(
    ncpus = 16,
    memory = 64,
    walltime = 36000,
    measure.memory = TRUE
)

main <- function() {
    if (dir.exists(batchtools)) {
        reg <- loadRegistry(batchtools, writeable = TRUE)
        removeRegistry(wait = 0)
    } 
    reg <- makeRegistry(batchtools, seed = 1)

    meta <- as_tibble(readRDS(in_seurat_meta))

    meta_subset <- meta %>% 
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-")) %>%
        filter((region == "preCG" & cluster_cell_type == "microglia"))
    
    subcluster_tbs <- meta_subset %>%
        group_by(region, cluster_cell_type) %>%
        group_nest(keep = TRUE)

    subcluster_wk <- subcluster_tbs %>%
        mutate(resample_cellids = map(data, function(data) {
            #resample_meta(data, 1, nsize = nrow(data))
            return(list(data$cell_ids))
        })) %>%
        unnest(resample_cellids)
    
    clearRegistry()
    ids <- batchMap(liger_cluster_sobj, 
        args = list(
            cell_ids = subcluster_wk$resample_cellids
        ), 
        more.args = list(
            filepath = in_seurat_rdata,
            liblist = liblist
        ))
    subcluster_wk$job.id <- getJobTable()$job.id
    submitJobs(ids, resources = RESOURCES)
    waitForJobs()
    saveRDS(subcluster_wk, file.path(out_path_base, "subcluster_resamples.rds"))
}
 
resample_meta <- function(meta, ntests = 5, nsize = floor(nrow(meta) * 0.8), cell_id_col = "cell_ids") {
    set.seed(27)
    return(map(1:ntests, function(i) {
        sample(meta[[cell_id_col]], nsize, replace = FALSE)
    }))
}

# Worker function: subset sobj to this instance's cell_ids, then convert to liger object and run clustering
liger_cluster_sobj <- function(cell_ids, filepath, liblist) {
    k <- 20
    lambda <- 0.5

    stopifnot("Duplicate cell ids cause errors in seuratToLiger -- don't use bootstrapping"=any(!duplicated(cell_ids)))
    l <- lapply(liblist, require, character.only = TRUE)

    envdata <- new.env()
    load(file.path(filepath), envdata) 
    #subset_sobj <- subset(envdata$liger_so, cells = cell_ids)
    subset_sobj <- subset(envdata$nd_so, cells = cell_ids)

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
