# lchen
# filter liger subclusters by max percent of library_id.

liblist <- c("tidyverse", "Seurat", "openxlsx")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

META_FILE <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../analysis/seurat_lchen/"
liger_meta <- readRDS(META_FILE)

meta <- liger_meta %>%
    as_tibble() %>%
    mutate(ct_subcluster = paste(cluster_cell_type, liger_clusters, sep = "-"))

meta <- meta %>%
    filter(cluster_cell_type == cell_type)

meta_libcluster <- meta %>%
    group_by(region, ct_subcluster, library_id) %>%
    summarize(clinical_dx = unique(clinical_dx), library_cluster_ct = n(), .groups = "drop") %>%
    group_by(region, ct_subcluster, clinical_dx) %>%
    mutate(cluster_total = sum(library_cluster_ct),
           library_cluster_pct = library_cluster_ct / cluster_total)
    
perdx_tb <- meta_libcluster %>%
    mutate(above_10umi = library_cluster_ct > 10) %>%
    group_by(region, ct_subcluster, clinical_dx) %>%
    summarize(libs_above_10umi = sum(above_10umi))

percluster_tb <- perdx_tb %>%
    group_by(region, ct_subcluster) %>%
    mutate(libs3_above_10umi = libs_above_10umi >= 3) %>%
    summarize(subcluster_3libs_above_10umi = sum(libs3_above_10umi) > 0)



saveRDS(meta_libcluster, file.path(OUT_DIR, "liger_subcluster_filtered_meta.rds"), compress = FALSE)
saveRDS(percluster_tb, file.path(OUT_DIR, "liger_subcluster_filtered_props.rds"))
write.xlsx(list(subcluster_ct = meta_libcluster, dx_summary = perdx_tb, cluster_summary = percluster_tb), file.path(OUT_DIR, "subcluster_percent.xlsx"))
