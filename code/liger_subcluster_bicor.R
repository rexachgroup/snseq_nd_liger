# Test nd score against library celltype counts per subcluster.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "WGCNA", "broom")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
in_seurat_liger <- "../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_bulk_meta <- "../resources/individual_nd.rds"
out_path_base <- "../analysis/seurat_lchen/liger_subcluster_bicor/"
clusters_exclude_file <- "../resources/subclusters_removed_byQC_final.xlsx"

main <- function() {
    dir.create(out_path_base, recursive = FALSE, showWarnings = FALSE)
    liger_meta_in <- readRDS(in_seurat_liger)
    bulk_meta <- readRDS(in_bulk_meta)
    excludes <- read_xlsx(clusters_exclude_file)
    
    nd_tb <- bulk_meta %>%
        mutate(Autopsy.ID = paste0("P", Autopsy.ID)) %>%
        group_by(Autopsy.ID, type) %>%
        slice_head(n = 1) %>%
        select(Autopsy.ID, type, score)
    
    library_celltype_counts_full <- liger_meta_in %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(library_id, ct_subcluster) %>%
        summarize(
            autopsy_id = unique(autopsy_id),
            clinical_dx = unique(clinical_dx),
            cell_type = unique(cell_type),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            cluster_ct = n()
        )

    dx_splits <- tibble(dx_filter = list(dx_only = c("AD", "PiD", "PSP"), all = c("AD", "PiD", "PSP", "Control")), dx_type = names(dx_filter))

    meta_nd <- inner_join(library_celltype_counts_full, nd_tb, by = c("autopsy_id" = "Autopsy.ID")) %>%
        filter(!is.na(type)) %>%
        group_by(ct_subcluster, type)

    dx_bicor_splits <- dx_splits %>%
        mutate(nd_bicor_tb = pmap(., function(...) {
                cr <- list(...)
                writeLines(str_glue("{paste(cr$dx_filter, collapse = ' ')}"))
                meta_nd_tb <- meta_nd %>%
                    filter(clinical_dx %in% cr$dx_filter) %>%
                    summarize(
                        bicor = bicor_trycatch(cluster_ct, score),
                        cell_type = unique(cell_type),
                        cor.test = cor_trycatch(cluster_ct, score),
                        cor.pval = cor.test$p.value,
                        .groups = "drop"
                    )
                meta_nd_tb <- meta_nd_tb %>%
                    group_by(cell_type) %>%
                    mutate(cor.fdr = p.adjust(cor.pval))
                return(meta_nd_tb)
            })
        )
    dx_bicor_pivot <- dx_bicor_splits %>%
        unnest(nd_bicor_tb) %>%
        select(ct_subcluster, dx_type, type, bicor, cor.pval, cor.fdr) %>%
        pivot_wider(id_cols = c("ct_subcluster"), names_from = c("dx_type", "type"), values_from = c("bicor", "cor.pval", "cor.fdr"), names_sep = "-")

    meta_test_vars <- c("age", "pmi")
    meta_bicor_splits <- dx_splits %>%
        mutate(nd_bicor_tb = pmap(., function(...) {
                cr <- list(...)
                meta_tb <- library_celltype_counts_full %>%
                    filter(clinical_dx %in% cr$dx_filter) %>%
                    group_by(ct_subcluster) %>%
                    summarize(
                        cell_type = unique(cell_type), 
                        across(all_of(meta_test_vars), function(x) { bicor_trycatch(x, cluster_ct) }, .names = "bicor_{.col}"), 
                        across(all_of(meta_test_vars), function(x) { cor_trycatch(x, cluster_ct)}, .names = "cor_{.col}"),
                        across(all_of(paste0("cor_", meta_test_vars)), function(x) { x$p.value }, .names = "pval_{.col}")
                    )
                meta_tb <- meta_tb %>%
                    group_by(cell_type) %>%
                    mutate(across(all_of(paste0("pval_cor_", meta_test_vars)), function(x) { p.adjust(x) }, .names = "fdr_{.col}"))
                return(meta_tb)
            }))

    meta_bicor_unnest <- meta_bicor_splits %>%
        unnest(nd_bicor_tb) %>%
        select(-starts_with("cor")) %>%
        select(ct_subcluster, dx_type, cell_type, all_of(contains(meta_test_vars))) %>%
        pivot_wider(id_cols = c("ct_subcluster"), names_from = c("dx_type"), values_from = contains(meta_test_vars), names_sep = "-")

    liger_meta <- liger_meta_in %>%
        group_by(ct_subcluster) %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster)
    nd_data <- mk_nd_data(liger_meta, nd_tb) %>%
        select(-cor.test)
    write_csv(dx_bicor_pivot, file.path(out_path_base, "ndscore_bicor_tests.csv"))
    write_csv(meta_bicor_unnest, file.path(out_path_base, "meta_bicor_tests.csv"))
    write_csv(nd_data, file.path(out_path_base, "ndscore_counts_bicor_tests.csv"))
}

mk_nd_data <- function(liger_meta, nd_tb) {
    library_celltype_counts_full <- liger_meta %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(library_id, ct_subcluster) %>%
        summarize(
            autopsy_id = unique(autopsy_id),
            clinical_dx = unique(clinical_dx),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            cluster_ct = n(),
        )

    meta_nd_tb <- left_join(library_celltype_counts_full, nd_tb, by = c("autopsy_id" = "Autopsy.ID")) %>%
        filter(!is.na(type)) %>%
        group_by(ct_subcluster, type) %>%
        summarize(
            bicor = bicor(cluster_ct, score, use = "pairwise.complete.obs", pearsonFallback = "none")[, 1],
            cor.test = tidy(cor.test(cluster_ct, score)),
            cor.pval = cor.test$p.value,
            cor.fdr = p.adjust(cor.pval),
            cor.signif = symnum_signif(cor.pval, cor.fdr)
        )

    return(meta_nd_tb)
}

symnum_signif <- function(pval, pval.fdr) {
    ifelse(
        pval.fdr < 0.1,
        symnum(pval.fdr, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", " ")),
        ifelse(
            pval < 0.05,
            "#",
            " "
        )
    )
}

bicor_trycatch <- function(x, y) { tryCatch(bicor(x, y, use = "pairwise.complete.obs", pearsonFallback = "none")[, 1], error = function(x) { return(tibble()) } ) }
cor_trycatch <- function(x, y) { tryCatch(tidy(cor.test(x, y)), error = function(x) { return (tibble()) }) }

if (!interactive()) main()
