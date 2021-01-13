liblist <- c("tidyverse", "batchtools")
l <- lapply(liblist, require, character.only = TRUE)

in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../analysis/pci_import/pci_seurat_meta.rds"
in_liger_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_clusters <- "../../analysis/seurat_lchen/cluster_stability_subjresample/reclusters.rds"
in_batchtools <- "../../analysis/seurat_lchen/cluster_stability_subjresample/batchtools/"
out_jaccard <- "../../analysis/seurat_lchen/cluster_stability_subjresample/recluster_jaccard.csv"
in_subcluster_wk <- "../../analysis/seurat_lchen/cluster_stability_subjresample/subcluster_resamples.rds"

main <- function() {
    meta <- readRDS(in_liger_meta)
    reg <- loadRegistry(in_batchtools)
    subcluster_wk <- readRDS(in_subcluster_wk)

    meta_subset <- meta %>% 
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-")) %>%
        filter((region == "preCG" & cluster_cell_type == "microglia") | (region == "insula" & cluster_cell_type == "excitatory"))

    subcluster_tbs <- meta_subset %>%
        group_by(region, cluster_cell_type) %>%
        group_nest(keep = TRUE)


    subcluster_tbs <- subcluster_tbs %>%
        mutate(job_id = pmap(., function(...) {
            current_row <- tibble(...)
            subcluster_wk %>% filter(region %in% current_row$region, cluster_cell_type %in% current_row$cluster_cell_type) %>% pluck("job.id")
        }))
        
    subcluster_tbs <- subcluster_tbs %>%
        mutate(liger_clustering = pmap(., function(job_id, ...) {
            map(job_id, extract_liger_reg)
        }))

    subcluster_tbs <- subcluster_tbs %>%
        mutate(jaccard_index = pmap(., function(data, liger_clustering, ...) {
            map(liger_clustering, function(x) {calc_liger_jacc_index(data, x)})
        }))


    jacc_unnest <- subcluster_tbs %>%
        select(-data, -liger_clustering) %>%
        unnest(c(jaccard_index, job_id)) %>%
        unnest(jaccard_index)
    
    jaccard_max <- jacc_unnest %>%
        group_by(region, cluster_cell_type, orig_cluster, job_id) %>%
        summarize(max_index_cluster = subs_cluster[which.max(jaccard_index)], max_index = max(jaccard_index), .groups = "drop") %>%
        arrange(job_id, orig_cluster)

    write_csv(jaccard_max, out_jaccard)
}

extract_liger_reg <- function(job.id, reg = getDefaultRegistry()) {
    lo <- loadResult(job.id, reg)
    return(tibble(cell_ids = rownames(lo@cell.data), liger_clusters = lo@clusters))
}


calc_liger_jacc_index <- function(orig, subsample) { 
    orig_cluster_ids <- unique(orig$liger_clusters)
    subs_cluster_ids <- unique(subsample$liger_clusters)
    cmp_tbl <- cross_df(list(orig_cluster = orig_cluster_ids, subs_cluster = subs_cluster_ids))
    cmp_tbl <- mutate(cmp_tbl, jaccard_index = pmap_dbl(cmp_tbl, function(...) {
        current_row <- tibble(...)
        orig_cell_ids <- orig %>% filter(liger_clusters == current_row$orig_cluster) %>% pluck("cell_ids")
        subs_cell_ids <- subsample %>% filter(liger_clusters == current_row$subs_cluster) %>% pluck("cell_ids")
        return(length(intersect(orig_cell_ids, subs_cell_ids)) / length(union(orig_cell_ids, subs_cell_ids)))
    }))
    return(cmp_tbl)
}

if (!interactive) {
    main()
}
