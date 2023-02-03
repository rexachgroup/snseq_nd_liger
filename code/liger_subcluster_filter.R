# Filter liger subclusters by minimum cell counts per library / dx.
# Also, filter by final subcluster filtering list to match plots.

liblist <- c("tidyverse", "Seurat", "openxlsx", "readxl")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

META_FILE <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
clusters_exclude_file <- "../resources/subclusters_removed_byQC_final.xlsx"
OUT_DIR <- "../analysis/seurat_lchen/"
liger_meta <- readRDS(META_FILE)
excludes <- read_xlsx(clusters_exclude_file)

meta <- liger_meta %>%
    as_tibble() %>%
    mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"))

meta <- meta %>%
    filter(cluster_cell_type == cell_type)

# Sum per library id + cluster.
meta_libcluster <- meta %>%
    group_by(region, ct_subcluster, library_id) %>%
    summarize(clinical_dx = unique(clinical_dx), library_cluster_ct = n(), .groups = "drop") %>%
    group_by(region, ct_subcluster, clinical_dx) %>%
    mutate(cluster_total = sum(library_cluster_ct),
           library_cluster_pct = library_cluster_ct / cluster_total)

# Flag and sum libraries above 10 cells per subcluster.
# Get passing library count per dx, then sum cells per dx.
perdx_tb <- meta_libcluster %>%
    mutate(above_10umi = library_cluster_ct > 10) %>%
    group_by(region, ct_subcluster, clinical_dx) %>%
    summarize(libs_above_10umi = sum(above_10umi), n = sum(library_cluster_ct))

# Count subclusters with more than 3 libraries per dx with > 10 cells.
# Propagate total cell count. 
percluster_tb <- perdx_tb %>%
    group_by(region, ct_subcluster) %>%
    mutate(libs3_above_10umi = libs_above_10umi >= 3) %>%
    summarize(subcluster_3libs_above_10umi = sum(libs3_above_10umi) > 0, n = sum(n))

# Drop subclusters by manual filter list to get matching list with plots.
perdx_tb_filt <- perdx_tb %>%
    filter(!ct_subcluster %in% excludes$ct_subcluster)
percluster_tb_filt <- percluster_tb %>%
    filter(!ct_subcluster %in% excludes$ct_subcluster)

saveRDS(meta_libcluster, file.path(OUT_DIR, "liger_subcluster_filtered_meta.rds"), compress = FALSE)
saveRDS(percluster_tb, file.path(OUT_DIR, "liger_subcluster_filtered_props.rds"))
write.xlsx(list(subcluster_ct = meta_libcluster, dx_summary = perdx_tb, cluster_summary = percluster_tb, dx_filt = perdx_tb_filt, cluster_filt = percluster_tb_filt), file.path(OUT_DIR, "subcluster_percent.xlsx"))
