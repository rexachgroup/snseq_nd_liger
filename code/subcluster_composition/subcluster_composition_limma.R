# limma across dx contrasts / dx splits within a celltype.
liblist <- c("tidyverse", "broom", "edgeR", "limma")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma/"
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
              liger_cluster = unique(liger_clusters),
              subcluster_ct = n())

library_celltype_counts_full <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(
        clinical_dx = unique(clinical_dx),
        age = unique(age),
        pmi = unique(pmi),
        sex = unique(sex),
        log_pmi = log(unique(pmi)),
        mean_percent_mito = mean(percent_mito),
        median_genes = median(number_genes),
        cluster_ct = n(),
    )

library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
    mutate(subcluster_pct_norm = subcluster_ct / cluster_ct)

library_cluster_mats <- library_pct %>%
    group_by(region, cell_type) %>%
    group_nest(keep = TRUE) %>%
    mutate(matrix = map(data, function(tb) { 
        tb %>%
            select(ct_subcluster, library_id, subcluster_ct) %>%
            pivot_wider(names_from = "library_id", values_from = "subcluster_ct", values_fill = 0) %>%
            column_to_rownames("ct_subcluster") %>%
            as.matrix
    }))



# for each data + matrix row, fit limma model using form. Compute contrasts + extract eBayes coefficients
limma_mutate <- function(tb, form) {
    tb %>%
    filter(map_lgl(matrix, ~nrow(.) > 2)) %>%
    mutate(limma = pmap(list(data, matrix), function(data, matrix){
        dgelist <- DGEList(matrix)
        dgelist_norm <- calcNormFactors(dgelist, method = "TMM")
        
        dx_fct <- data %>% group_by(library_id) %>% 
            summarize(
                dx = make.names(unique(clinical_dx)),
                pmi = unique(pmi),
                log_pmi = unique(log_pmi),
                sex = unique(sex),
                age = unique(age),
                mean_percent_mito = unique(mean_percent_mito),
                median_genes = unique(median_genes),
                .groups = "drop"
            )
        design  <- model.matrix(as.formula(form), data = dx_fct)
        
        elist_voom <- voom(dgelist_norm, design)
        dx_fit <- lmFit(elist_voom, design)

        limma_contrasts <- makeContrasts(
            contrasts = c("dxAD - dxControl",
            "dxbvFTD - dxControl",
            "dxPSP.S - dxControl",
            "(dxAD + dxbvFTD + dxPSP.S)/3 - dxControl",
            "dxAD - (dxControl + dxbvFTD + dxPSP.S) / 3",
            "dxbvFTD - (dxControl + dxAD + dxPSP.S) / 3",
            "dxPSP.S - (dxControl + dxAD + dxbvFTD) / 3"
            ),
            levels = colnames(coef(dx_fit))
        )
        dx_contrast_fit <- contrasts.fit(dx_fit, limma_contrasts)
        return(dx_contrast_fit)
    })) %>%
    mutate(
        limma_pval = map(limma, function(limma) {
            tryCatch({
                dx_bayes <- eBayes(limma)
                dx_bayes$p.value %>%
                    as.data.frame %>%
                    rename_with(function(x) paste0("pval", x), everything()) %>%
                    as_tibble(rownames = "ct_subcluster")
            }, error = function(x) NA)
        }),
        limma_beta = map(limma, function(limma) {
            limma$coefficients %>%
                as.data.frame %>%
                rename_with(function(x) paste0("beta", x), everything()) %>%
                as_tibble(rownames = "ct_subcluster")
        }),
        formula = form
    )
}

dx_splits <- list(
    ad_v_ctl = list("AD", "Control"),
    bvftd_v_ctl = list("bvFTD", "Control"),
    psps_v_ctl = list("PSP-S", "Control"),
    alldx_v_ctl = list(c("AD", "bvFTD", "PSP-S"), "Control"),
    ad_v_all = list("AD", c("bvFTD", "PSP-S", "Control")),
    bvftd_v_all = list("bvFTD", c("AD", "PSP-S", "Control")),
    psps_v_all = list("PSP-S", c("AD", "bvFTD", "Control"))
)

# generate_group_splits <- functi(factor_name, factor_split_list, tb) {
#     imap(factor_split_list, function(split, split_name) {
#         group1 <- split[[1]]
#         group2 <- split[[2]]
#         group1_rows <- tb[[factor_name]] %in% group1
#         group2_rows <- tb[[factor_name]] %in% group2
#         subset_tb <- tb[group1_rows | group2_rows, ]
# 
#         subset_tb <- subset_tb %>%
#             mutate(
#                 group = split_name,
#                 group_split = ifelse(.data[[factor_name]] %in% group1, 1, ifelse(.data[[factor_name]] %in% group2, 0, NA))
#             )
#         return(subset_tb)
#     })
# }
# 
# limma_mutate_split <- function(tb, form) { 
#     tb <- tb %>%
#     filter(map_lgl(matrix, ~nrow(.) > 2)) %>%
#     mutate(limma = pmap(list(data, matrix), function(data, matrix){
#         data_subsets <- generate_group_splits("clinical_dx", dx_splits, data)
# 
#         data_subsets[lapply(data_subsets, nrow) > 0] %>%
#             map(function(data_subset) {
#                 dxtb <- data_subset %>% group_by(library_id) %>%
#                     summarize(
#                         dx = make.names(unique(clinical_dx)),
#                         group_split = unique(group_split),
#                         pmi = unique(pmi),
#                         log_pmi = unique(log_pmi),
#                         sex = unique(sex),
#                         age = unique(age),
#                         mean_percent_mito = unique(mean_percent_mito),
#                         median_genes = unique(median_genes),
#                         .groups = "drop"
#                     )
#                 dxmatrix <- matrix[, dxtb$library_id]
#                 
#                 dx_fit <- NA
#                 tryCatch({
#                     dgelist <- DGEList(dxmatrix)
#                     dgelist_norm <- calcNormFactors(dgelist, method = "TMM")
#                     design  <- model.matrix(as.formula(form), data = dxtb) 
#                     elist_voom <- voom(dgelist_norm, design)
#                     dx_fit <- suppressWarnings(lmFit(elist_voom, design))
#                 }, 
#                 error = function(x) {},
#                 warning = function(x) {},
#                 finally = return(dx_fit))
#             })
#     }))
#     tb <- tb %>%
#     mutate(
#         limma_pval = map(limma, function(limma_list) {
#             tryCatch({
#             pval_list <- map(limma_list, function(lim) {
#                 dx_bayes <- eBayes(lim)
#                 dx_bayes$p.value %>%
#                     as_tibble(rownames = "ct_subcluster") %>%
#                     rename(group_split_pval = group_split)
#             })
# 
#             pval_list[!is.na(pval_list)] %>%
#                 bind_rows(.id = "group_split_name")
#             }, error = function(x) {})
#         }),
#         formula = form
#     )
# }


summarize_limma <- function(tb) {
    tb %>%
        filter(!is.na(limma_pval)) %>%
        mutate(limma_cmb = pmap(list(limma_beta, limma_pval), ~inner_join(.x, .y, by = "ct_subcluster"))) %>%
        unnest(limma_cmb) %>%
        glimpse
}

write_limma <- function(tb, path) {
    tb %>% 
        filter(!is.na(limma_pval)) %>%
        mutate(limma_cmb = pmap(list(limma_beta, limma_pval), ~inner_join(.x, .y, by = "ct_subcluster"))) %>%
        unnest(limma_cmb) %>%
        select(-data, -matrix, -limma, -limma_pval, -limma_beta) %>%
        write_tsv(path)
}


#limma_dx <- limma_mutate(library_cluster_mats, form = "~0 + dx")
#limma_dx_log_pmi <- limma_mutate(library_cluster_mats, form = "~0 + dx + log_pmi + age + sex + mean_percent_mito + median_genes")
limma_dx_pmi <- limma_mutate(library_cluster_mats, form = "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes")
#limma_dx_split <- limma_mutate_split(library_cluster_mats, form = "~0 + group_split")

# summarize_limma(limma_dx)
# summarize_limma(limma_dx_log_pmi)
summarize_limma(limma_dx_pmi)

# write_limma(limma_dx, file.path(OUT_DIR, "subcluster_composition_dx.tsv"))
# write_limma(limma_dx_log_pmi, file.path(OUT_DIR, "subcluster_composition_dx_log_pmi.tsv"))
write_limma(limma_dx_pmi, file.path(OUT_DIR, "subcluster_composition_dx_pmi.tsv"))
#write_limma(limma_dx_split, file.path(OUT_DIR, "subcluster_composition_dx_split.tsv"))

# saveRDS(limma_dx, file.path(OUT_DIR, "subcluster_composition_dx.rds"))
# saveRDS(limma_dx_log_pmi, file.path(OUT_DIR, "subcluster_composition_dx_log_pmi.rds"))
saveRDS(limma_dx_pmi, file.path(OUT_DIR, "subcluster_composition_dx_pmi.rds"))

