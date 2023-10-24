# Test nd score against library celltype counts per subcluster.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "WGCNA", "broom")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
SUBCLUSTER_FILTERED <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/01_subcluster_filtered_meta.csv"
in_nd_tb <- "../../analysis/bulk_meta/nd_tb.rds"
out_path_base <- "../../analysis/seurat_lchen/subcluster_excluded_dge/celltype-nd-cor/"

main <- function() {
    dir.create(out_path_base, recursive = FALSE, showWarnings = FALSE)
    liger_meta_in <- read_csv(SUBCLUSTER_FILTERED)
    nd_tb <- readRDS(in_nd_tb)

    library_celltype_counts_full <- liger_meta_in %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(library_id, ct_subcluster) %>%
        summarize(
            autopsy_id = unique(autopsy_id),
            clinical_dx = unique(clinical_dx),
            region = unique(region),
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

    liger_meta_dxonly <- liger_meta_in %>%
        group_by(ct_subcluster) %>%
        filter(clinical_dx != "Control")

    nd_data_ct <- liger_meta_dxonly %>%
        group_by(cluster_cell_type) %>%
        group_nest(keep = TRUE) %>%
        mutate(nd_data = map(data, mk_nd_data, nd_tb))
    nd_data_ct_region <- liger_meta_dxonly %>%
        group_by(cluster_cell_type, region) %>%
        group_nest(keep = TRUE) %>%
        mutate(nd_data = map(data, mk_nd_data, nd_tb))

    nd_data_ct_region_dx <- liger_meta_dxonly %>%
        group_by(cluster_cell_type, region, clinical_dx) %>%
        group_nest(keep = T) %>%
        mutate(nd_data = map(data, mk_nd_data, nd_tb))

    nd_data_ct %>%
        unnest(nd_data) %>%
        select(-data, -cor.test) %>%
        write_csv(file.path(out_path_base, "fdr_celltype_ndscore_bicor_tests_dxonly.csv"))

    nd_data_ct_region %>%
        unnest(nd_data) %>%
        select(-data, -cor.test) %>%
        write_csv(file.path(out_path_base, "fdr_celltype_region_ndscore_counts_bicor_tests_dxonly.csv"))
    
    nd_data_ct_region_dx %>%
        unnest(nd_data) %>%
        select(-data, -cor.test) %>%
        write_csv(file.path(out_path_base, "fdr_celltype_dx_region_ndscore_counts_bicor_tests.csv"))

    write_csv(meta_bicor_unnest, file.path(out_path_base, "meta_bicor_tests.csv"))
}

mk_nd_data <- function(liger_meta, nd_tb) {
    library_celltype_counts_full <- liger_meta %>%
        mutate(ct_subcluster = fct_drop(ct_subcluster)) %>%
        group_by(library_id, ct_subcluster) %>%
        summarize(
            autopsy_id = unique(autopsy_id),
            clinical_dx = unique(clinical_dx),
            region = unique(region),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            cluster_ct = n(),
            .groups = "drop"
        )

    # bicor / cor.test on cell counts per subcluster.
    meta_nd_tb <- inner_join(library_celltype_counts_full, nd_tb, by = c("library_id"), multiple = "all", relationship = "many-to-many") %>%
        filter(!is.na(type)) %>%
        group_by(ct_subcluster, type) %>%
        summarize(
            bicor = bicor_trycatch(cluster_ct, score),
            cor.test = cor_trycatch(cluster_ct, score),
            .groups = "drop"
        )

    if (nrow(meta_nd_tb) == 0) { return(tibble()) }
    # do fdr correction per nd type.
    meta_nd_tb <- meta_nd_tb %>%
        group_by(type) %>%
        mutate(
            cor.pval = cor.test$p.value,
            cor.fdr = p.adjust(cor.pval),
            cor.n = n(),
            cor.signif = symnum_signif(cor.pval, cor.fdr)
        )
    # fill in subclusters that didn't have a sample in nd_tb
    meta_nd_tb <- meta_nd_tb %>% ungroup %>% complete(expand(meta_nd_tb, ct_subcluster, type))

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
