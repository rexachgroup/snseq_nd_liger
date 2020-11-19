# Lawrence Chen 2020-11-03
# anova on cell proportion changes in all data.

liblist <- c("multcomp", "tidyverse", "broom", "patchwork")
lapply(liblist, require, character.only = TRUE)

#in_seurat_rds <- "../analysis/pci_import/20201028/tables/pci_seurat.rds"
META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/anova/"
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
    summarize(cell_type = unique(cell_type), clinical_dx = unique(clinical_dx), subcluster_ct = n()) %>%
    print(n = 50)

dx_splits <- list(
    ad_v_ctl = list("AD", "Control"),
    bvftd_v_ctl = list("bvFTD", "Control"),
    psps_v_ctl = list("PSP-S", "Control"),
    alldx_v_ctl = list(c("AD", "bvFTD", "PSP-S"), "Control"),
    ad_v_all = list("AD", c("bvFTD", "PSP-S", "Control")),
    bvftd_v_all = list("bvFTD", c("AD", "PSP-S", "Control")),
    psps_v_all = list("PSP-S", c("AD", "bvFTD", "Control"))
)

dx_contrast_mat <- cbind(
    ad_v_ctl = c(-1, 1, 0, 0), # ad vs. ctl
    bvftd_v_ctl = c(-1, 0, 1, 0), # bvftd vs. ctl
    psps_v_ctl = c(-1, 0, 0, 1), # psp-s vs. ctl
    alldx_v_ctl = c(-1, 1/3, 1/3, 1/3), # all dx vs. ctl
    ad_v_all = c(-1/3, 1, -1/3, -1/3), # ad vs. all
    bvftd_v_all = c(-1/3, -1/3, 1, -1/3), # bvftd vs. all
    psps_v_all = c(-1/3, -1/3, -1/3, 1) # psp-s vs. all 
)

#contrast_obj <- mcp(clinical_dx = dx_contrast_mat[, 1])
contrast_obj <- mcp(clinical_dx = c(
    ad_v_ctl = "AD - Control == 0",
    bvftd_v_ctl = "bvFTD - Control == 0",
    psps_v_ctl = "`PSP-S` - Control == 0",
    alldx_v_ctl = "(AD + bvFTD + `PSP-S`)/3 - Control == 0",
    ad_v_all = "AD - (bvFTD + `PSP-S` + Control)/3 == 0",
    bvftd_v_all = "bvFTD - (AD + `PSP-S` + Control)/3 == 0",
    psps_v_all = "`PSP-S` - (AD + bvFTD + Control)/3 == 0"
))

# because I can't figure out how contrasts work, we're subsetting by each of the dx splits.
generate_group_splits <- function(factor_name, factor_split_list, tb) {
    imap(factor_split_list, function(split, split_name) {
        group1 <- split[[1]]
        group2 <- split[[2]]
        group1_rows <- tb[[factor_name]] %in% group1
        group2_rows <- tb[[factor_name]] %in% group2
        subset_tb <- tb[group1_rows | group2_rows, ]

        subset_tb <- subset_tb %>%
            mutate(
                group = split_name,
                group_split = ifelse(.data[[factor_name]] %in% group1, 1, ifelse(.data[[factor_name]] %in% group2, 0, NA))
            )
        return(subset_tb)
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
                lm_list %>%
                    map(anova) %>%
                    map(tidy) %>%
                    bind_rows(.id = "group") %>%
                    pivot_wider(
                        names_from = "term",
                        values_from = c("df", "sumsq", "meansq", "statistic", "p.value")
                    ) %>%
                    glimpse
        }),
        aov_form = {{aov_form}}
    ) %>%
    unnest(aov_tidy)
}

# 1. Count percentages cells for the seurat cell type

library_celltype_counts_full <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(cluster_ct = n()) %>%
    print(n = 50)

library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct) %>%
    glimpse()

aov_form <- "subcluster_pct_norm ~ group_split"
library_cell_pct_anova <- aov_dx_subsets(library_pct, aov_form)
    

# 2. Count percentages for the seurat cell type after filtering subclusters

library_celltype_counts_filter <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    filter(ct_subcluster %in% percluster_tb$ct_subcluster) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(cluster_ct = n()) %>%
    print(n = 50)

library_pct_filter <- inner_join(library_subcluster_counts, library_celltype_counts_filter, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct) %>%
    glimpse()

library_cell_pct_anova_filter <- aov_dx_subsets(library_pct_filter, aov_form)


# 3. glht contrast testing

library_pct_glht <- library_pct_filter %>% group_by(ct_subcluster) %>%
    group_nest() %>%
    mutate(
        aov_objs = map(data, ~aov(as.formula("subcluster_pct_norm ~ clinical_dx"), .)),
        glht_objs = map(aov_objs, ~glht(., linfct = contrast_obj)),
        glht_tidy = map(glht_objs, tidy)
    )

# 4. data with subcluster counts 
library_cell_ct_anova <- aov_dx_subsets(library_pct_filter, "subcluster_ct ~ group_split")

# 5. ggplot count distribution
library_pct_patchwork <- library_pct_filter %>%
    group_by(ct_subcluster) %>%
    group_split() %>%
    map(~ggplot(data = ., aes(x = clinical_dx, y = subcluster_pct_norm, color = clinical_dx)) + geom_boxplot() + ggtitle(unique(.$ct_subcluster)))
    
library_ct_patchwork <- library_pct_filter %>%
    group_by(ct_subcluster) %>%
    group_split() %>%
    map(~ggplot(data = ., aes(x = clinical_dx, y = subcluster_ct, color = clinical_dx)) + geom_boxplot() + ggtitle(unique(.$ct_subcluster)))

# Output
write_csv(library_celltype_counts_full, file.path(OUT_DIR, "library_celltype_counts_full.csv"))
write_csv(library_celltype_counts_filter, file.path(OUT_DIR, "library_celltype_counts_filter.csv"))

saveRDS(library_cell_pct_anova, file.path(OUT_DIR, "library_cell_pct_anova.rds"))
library_cell_pct_anova %>%
    dplyr::select(-data, -dx_subsets) %>%
    write_csv(file.path(OUT_DIR, "library_cell_pct_anova.csv"))

saveRDS(library_cell_pct_anova_filter, file.path(OUT_DIR, "library_cell_pct_anova_filter.rds"))
library_cell_pct_anova_filter %>%
    dplyr::select(-data, -dx_subsets) %>%
    write_csv(file.path(OUT_DIR, "library_cell_pct_anova_filter.csv"))

saveRDS(library_pct_glht, file.path(OUT_DIR, "library_pct_glht.rds"))
library_pct_glht %>%
    dplyr::select(ct_subcluster, glht_tidy) %>%
    unnest(glht_tidy) %>%
    write_csv(file.path(OUT_DIR, "library_pct_glht.csv"))

saveRDS(library_cell_ct_anova, file.path(OUT_DIR, "library_cell_ct_anova.rds"))
library_cell_ct_anova %>%
    dplyr::select(-data, -dx_subsets) %>%
    write_csv(file.path(OUT_DIR, "library_cell_ct_anova.csv"))

pdf(file.path(OUT_DIR, "library_cell_pct_boxplot.pdf"), width = 5, height = 5)
tryCatch(print(library_pct_patchwork))
dev.off()

pdf(file.path(OUT_DIR, "library_cell_ct_boxplot.pdf"), width = 5, height = 5)
tryCatch(print(library_ct_patchwork))
dev.off()
