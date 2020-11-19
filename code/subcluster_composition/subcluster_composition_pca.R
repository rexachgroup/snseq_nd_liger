## PCA of library subcluster counts.
liblist <- c("multcomp", "tidyverse", "broom", "patchwork", "GGally")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/pca/"
SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

meta <- readRDS(META_FILE)
percluster_tb <- readRDS(SUBCLUSTER_FILTER_FILE) %>% 
    filter(subcluster_3libs_above_10umi) %>%
    mutate(ct_subcluster = paste(region, ct_subcluster, sep = "-"))

meta <- meta %>%
    mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
        clinical_dx = fct_relevel(clinical_dx, "Control"))

ctl_filter <- meta %>% group_by(ct_subcluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)

library_subcluster_counts <- meta %>%
    filter(cluster_cell_type == cell_type,
           ct_subcluster %in% percluster_tb$ct_subcluster,
           ct_subcluster %in% ctl_filter$ct_subcluster) %>%
    group_by(region, library_id, ct_subcluster) %>%
    summarize(cell_type = unique(cell_type), 
              clinical_dx = unique(clinical_dx),
              age = unique(age),
              pmi = unique(pmi),
              percent_mito = mean(percent_mito),
              median_genes = median(number_genes),
              mean_reads = mean(number_umi),
              log_pmi = unique(log(pmi)),
              log_mito = unique(log(percent_mito)),
              subcluster_ct = n())

library_celltype_counts_full <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(cluster_ct = n())

library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct)

library_cluster_mats <- library_pct %>%
    group_by(region, cell_type) %>%
    group_nest(keep = TRUE) %>%
    mutate(matrix = map(data, function(tb) { 
        tb %>%
            select(ct_subcluster, library_id, subcluster_pct_norm) %>%
            pivot_wider(names_from = "library_id", values_from = "subcluster_pct_norm", values_fill = 0) %>%
            column_to_rownames("ct_subcluster") %>%
            as.matrix
    }))

library_cluster_mats <- library_cluster_mats %>%
    mutate(prcmp = map(matrix, function(x) x %>% t %>% prcomp))

library_cluster_mats <- library_cluster_mats %>%
    mutate(cluster_gg = pmap(list(data, prcmp), function(data, prcmp) {
        pcs <- prcmp$rotation[, 1:ifelse(nrow(prcmp$rotation) < 5, nrow(prcmp$rotation), 5)] %>%
            as_tibble(rownames = "ct_subcluster")
        ggdat <- data %>% left_join(pcs, by = "ct_subcluster") %>%
            dplyr::select(-"library_id", -"ct_subcluster")
        glimpse(ggdat)
        ggpairs(ggdat, aes(color = clinical_dx), lower = list(combo = wrap("facethist", bins = 30))) + ggtitle(paste0(unique(data$region), unique(data$cell_type)))
    }))

pdf(file.path(OUT_DIR, "pca.pdf"), width = 20, height = 20)
walk(library_cluster_mats$cluster_gg, ~tryCatch(print(.), warning = function(x) print(x), error = function(x) print(x)))
dev.off()

