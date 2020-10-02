# lchen
# filter liger subclusters by max percent of library_id.

liblist <- c("tidyverse", "Seurat")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

SEURAT_FILE <- "../analysis/seurat_lchen/liger_subcluster_merged.rds"
OUT_DIR <- "../analysis/seurat_lchen/"

liger_so <- readRDS(SEURAT_FILE)
write_csv(liger_so@meta.data, file.path(OUT_DIR, "liger_subcluster_meta.csv"))

meta <- liger_so@meta.data
liger_so <- AddMetaData(liger_so,
    paste(meta[["cluster_cell_type"]], meta[["liger_clusters"]], sep = "-"),
    col.name = "ct_subcluster")

meta <- liger_so@meta.data %>%
    as_tibble(rownames = "UMI")

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

liger_subset <- subset(liger_so, cells = meta %>% filter(ct_subcluster %in% subcluster_filtered) %>% pull(UMI))
saveRDS(OUT_DIR, file.path(OUT_DIR, "liger_subcluster_subset.rds"), compress = FALSE)
