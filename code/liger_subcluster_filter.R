# lchen
# filter liger subclusters by max percent of library_id.

liblist <- c("tidyverse", "Seurat")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

META_FILE <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../analysis/seurat_lchen/"
liger_meta <- readRDS(META_FILE)

meta <- liger_meta %>%
    mutate(ct_subcluster = paste(cluster_cell_type, liger_clusters, sep = "-"))

subcluster_library_props <- meta %>%
    group_by(ct_subcluster, library_id) %>%
    summarize(umi_ct = n(), .groups = "drop") %>%
    group_by(ct_subcluster) %>%
    mutate(pct_of_subcluster = umi_ct / sum(umi_ct))

write_csv(subcluster_library_props, file.path(OUT_DIR, "liger_subcluster_counts.csv"))

subcluster_filtered <- subcluster_library_props %>%
    group_by(ct_subcluster) %>%
    summarize(max_lib_subcluster_pct = max(pct_of_subcluster), .groups = "drop") %>%
    filter(max_lib_subcluster_pct < 0.5) %>%
    pull(ct_subcluster)

meta <- meta %>%
    filter(ct_subcluster %in% subcluster_filtered)

saveRDS(meta, file.path(OUT_DIR, "liger_subcluster_subset.rds"), compress = FALSE)
