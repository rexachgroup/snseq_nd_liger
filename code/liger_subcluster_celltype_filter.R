# Filter based on ratio of majority cell type in each subcluster.
set.seed(0)
liblist <- c("Seurat", "tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <-
  "../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../analysis/seurat_lchen/liger_subcluster_metadata.rds"

## Outputs
out_path_base <- "../analysis/seurat_lchen/liger_subcluster_celltype_filter/"

main <- function() {
    liger_meta <- readRDS(in_seurat_liger) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_subcluster = factor(ct_subcluster, levels = str_sort(unique(ct_subcluster), numeric = TRUE))
        )
 
    meta_subset <- liger_meta %>% 
        group_by(ct_subcluster)
    
    meta_subset_celltype_match <- liger_meta %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(ct_subcluster)
    
    meta_sum <- meta_subset %>%
        group_by(ct_subcluster) %>%
        summarize(n_before = n())

    meta_sum_celltype_match <- meta_subset_celltype_match %>%
        group_by(ct_subcluster) %>%
        summarize(n_after = n())
    
    meta_sum_ratio <- inner_join(meta_sum, meta_sum_celltype_match, by = "ct_subcluster") %>%
        mutate(ratio = n_after / n_before)

    subcluster_filtered <- meta_sum_ratio %>%
        filter(ratio > 0.7) %>%
        pluck("ct_subcluster")

    subcluster_override <- liger_meta %>%
        filter(region == "preCG" & cluster_cell_type == "inhibitory") %>%
        pluck("ct_subcluster") %>%
        unique

    subcluster_vec <- union(subcluster_filtered, subcluster_override)
    saveRDS(subcluster_vec, file.path(out_path_base, "celltype_filter.rds"))
}

if (!interactive()) {
    main()
}
