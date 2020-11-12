# Lawrence Chen 2020-11-03
# anova on cell proportion changes in all data.

liblist <- c("tidyverse", "broom")
lapply(liblist, require, character.only = TRUE)

#in_seurat_rds <- "../analysis/pci_import/20201028/tables/pci_seurat.rds"
META_FILE <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../analysis/clinical_dx_anova/"
SUBCLUSTER_FILTER_FILE <- "../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
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
    summarize(cell_type = unique(cell_type), clinical_dx = unique(clinical_dx), subcluster_ct = n()) %>%
    print(n = 50)

pairwise_subsets <- function(factor_name, tb) {
    baselvl <- levels(tb[[factor_name]])[1]
    subsets <- lapply(levels(tb[[factor_name]])[-1], function(testlvl) {
        intest <- (tb[[factor_name]] == baselvl) | (tb[[factor_name]] == testlvl)
        tb_subset <- tb[intest, ]
        tb_subset[[factor_name]] <- droplevels(tb_subset[[factor_name]]) 
        tb_subset[["subset"]] <- paste0(testlvl, factor_name)
        return(tb_subset)
    })
    return(setNames(subsets, nm = levels(tb[[factor_name]])[-1]))
}

dx_splits <- list(
    ad_v_ctl = c("AD", "Control"),
    bvftd_v_ctl = c("bvFTD", "Control"),
    psps_v_ctl = c("PSP-S", "Control"),
    alldx_v_ctl = c(c("AD", "bvFTD", "PSP-S"), "Control"),
    ad_v_all = c("AD", c("bvFTD", "PSP-S", "Control")),
    bvftd_v_all = c("bvFTD", c("AD", "PSP-S", "Control")),
    psps_v_all = c("PSP-S", c("AD", "bvFTD", "Control"))
)

dx_contrast_list <- cbind(
    ad_v_ctl = c(-1, 1, 0, 0), # ad vs. ctl
    bvftd_v_ctl = c(-1, 0, 1, 0), # bvftd vs. ctl
    psps_v_ctl = c(-1, 0, 0, 1), # psp-s vs. ctl
    alldx_v_ctl = c(-1, 1/3, 1/3, 1/3), # all dx vs. ctl
    ad_v_all = c(-1/3, 1, -1/3, -1/3), # ad vs. all
    bvftd_v_all = c(-1/3, -1/3, 1, -1/3), # bvftd vs. all
    psps_v_all = c(-1/3, -1/3, -1/3, 1) # psp-s vs. all 
)

# because I can't figure out how contrasts work, we're subsetting by each of the dx splits.
generate_group_splits <- function(factor_name, factor_split_list, tb) {
    imap(factor_split_list, function(split, split_name) {
        group1 <- split[[1]]
        group2 <- split[[2]]
        group1_rows <- tb[[factor_name]] %in% group1
        group2_rows <- tb[[factor_name]] %in% group2
        tb[["group"]] <- split_name
        tb$group_split <- vector(length = nrow(tb))
        tb$group_split[group1_rows] <- 1
        tb$group_split[group2_rows] <- 0
        return(tb[group1_rows | group2_rows, ])
    })
}

aov_dx_subsets <- function(tb, aov_form) {
    tb %>% 
    group_by(ct_subcluster) %>%
    group_nest(keep = TRUE) %>%
    mutate(
        dx_subsets = map(data, function(tb) {
            aov_dx_subsets <- generate_group_splits("clinical_dx", dx_splits, tb)
            aov_dx_subsets[lapply(aov_dx_subsets, nrow) > 0] %>%
                map(function(dxtb) {
                    lm(as.formula(aov_form), dxtb)
                })
        }),
        aov_tidy = map(dx_subsets, function(lm_list) {
            valid_lm <- map_lgl(lm_list, ~!is.null(.))
            if (any(valid_lm))
                lm_list[valid_lm] %>%
                    map(aov) %>%
                    map(tidy) %>%
                    bind_rows(.id = "group") %>%
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

aov_form <- "subcluster_pct_norm ~ group_split"
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
