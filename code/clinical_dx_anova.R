# Lawrence Chen 2020-11-03
# anova on cell proportion changes in all data.

liblist <- c("tidyverse", "broom")
lapply(liblist, require, character.only = TRUE)

#in_seurat_rds <- "../analysis/pci_import/20201028/tables/pci_seurat.rds"
META_FILE <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../analysis/clinical_dx_anova/"
SUBCLUSTER_FILTER_FILE <- "../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

#nd_so <- readRDS(in_seurat_rds)

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
    summarize(cell_type = unique(cell_type), clinical_dx = unique(clinical_dx), subcluster_ct = n()) %>%
    print(n = 50)

pairwise_subsets <- function(factor_name, tb) {
    baselvl <- levels(tb[[factor_name]])[1]
    subsets <- lapply(levels(tb[[factor_name]])[-1], function(testlvl) {
        intest <- tb[[factor_name]] == baselvl | tb[[factor_name]] == testlvl
        tb_subset <- tb[intest, ]
        tb_subset[[factor_name]] <- droplevels(tb_subset[[factor_name]]) 
        tb_subset[["subset"]] <- paste0(testlvl, factor_name)
        return(tb_subset)
    })
    return(setNames(subsets, nm = levels(tb[[factor_name]])[-1]))
}

aov_dx_subsets <- function(tb, aov_form) {
    tb %>% 
    group_by(ct_subcluster) %>%
    group_nest(keep = TRUE) %>%
    mutate(
        dx_subsets = map(data, function(tb) {
            aov_dx_subsets <- pairwise_subsets("clinical_dx", tb)
            aov_dx_subsets[lapply(aov_dx_subsets, nrow) > 0] %>%
                map(function(dxtb) {
                    if (length(levels(dxtb$clinical_dx)) == 2) {
                        aov(as.formula(aov_form), dxtb)
                    }
                })
        }),
        aov_tidy = map(dx_subsets, function(aov_list) {
            valid_aov <- map_lgl(aov_list, ~!is.null(.))
            if (any(valid_aov))
                aov_list[valid_aov] %>%
                    map(tidy) %>%
                    bind_rows(.id = "dx") %>%
                    pivot_wider(
                        names_from = "term",
                        values_from = c("df", "sumsq", "meansq", "statistic", "p.value")
                    ) %>%
                    glimpse
        })
    ) %>%
    unnest(aov_tidy)
}

# 1. Count all cells for the seurat cell type

library_celltype_counts_full <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(cluster_ct = n()) %>%
    print(n = 50)

library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct) %>%
    glimpse()

aov_form <- "subcluster_pct_norm ~ clinical_dx"
library_cell_counts_anova <- aov_dx_subsets(library_pct, aov_form)
    

# 2. Count cells for the seurat cell type after filtering subclusters

library_celltype_counts_filter <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    filter(ct_subcluster %in% percluster_tb$ct_subcluster) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(cluster_ct = n()) %>%
    print(n = 50)

library_pct_filter <- inner_join(library_subcluster_counts, library_celltype_counts_filter, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct) %>%
    glimpse()

aov_form <- "subcluster_pct_norm ~ clinical_dx"
library_cell_counts_anova_filter <- aov_dx_subsets(library_pct_filter, aov_form)

# 3. Output
write_csv(library_celltype_counts_full, file.path(OUT_DIR, "library_celltype_counts_full.csv"))
write_csv(library_celltype_counts_filter, file.path(OUT_DIR, "library_celltype_counts_filter.csv"))

saveRDS(library_cell_counts_anova, file.path(OUT_DIR, "library_cell_counts_anova.rds"))
library_cell_counts_anova %>%
    select(-data, -dx_subsets) %>%
    write_csv(file.path(OUT_DIR, "library_cell_counts_anova.csv"))

saveRDS(library_cell_counts_anova_filter, file.path(OUT_DIR, "library_cell_counts_anova_filter.rds"))
library_cell_counts_anova_filter %>%
    select(-data, -dx_subsets) %>%
    write_csv(file.path(OUT_DIR, "library_cell_counts_anova_filter.csv"))
