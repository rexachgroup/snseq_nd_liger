# Test contrasts on all limma metadata.
liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl")
lapply(liblist, require, character.only = TRUE)
set.seed(NULL)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma_metadata/"

main <- function() {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    in_meta <- readRDS(META_FILE)

    meta <- in_meta

    library_meta <- meta %>%
        group_by(library_id) %>%
        summarize(
            dx = make.names(unique(clinical_dx)),
            region = unique(region),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes)
        )

    ctl_filter <- meta %>% group_by(ct_subcluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)
    
    library_subcluster_counts <- meta %>%
        filter(
            cluster_cell_type == cell_type,
            ct_subcluster %in% ctl_filter$ct_subcluster
        ) %>%
        group_by(region, library_id, ct_subcluster) %>%
        summarize(
            cell_type = unique(cell_type), 
            liger_cluster = unique(liger_clusters),
            subcluster_ct = n()
        ) %>%
        group_by(region, cell_type) %>%
        group_nest(keep = TRUE)

    library_subcluster_meta_limma <- library_subcluster_counts %>%
        mutate(limma_meta_results = map(data, function(.x) {
            library_subset_meta <- filter(library_meta, library_id %in% .x$library_id)
            fit <- limma_dx_meta_fit(
                library_subset_meta,
                .x,
                library_name = "library_id",
                count_name = "subcluster_ct",
                class_name = "ct_subcluster",
                formula = "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes",
                randomize = FALSE
            )
            return(fit)
        }))


    library_subcluster_meta_tb <- library_subcluster_meta_limma %>%
        mutate(coefs = map(limma_meta_results, limma_extract_coefs)) %>%
        select(-data, -limma_meta_results) %>%
        unnest(coefs)
    write_csv(library_subcluster_meta_tb, file.path(OUT_DIR, "subcluster_composition_meta.csv"))
}

limma_extract_coefs <- function(limma_result) {
    limma_pval <- limma_result$bayes$p.value %>%
        as.data.frame
    limma_beta <- limma_result$fit$coefficients %>%
        as.data.frame
    limma_sigma <- limma_result$bayes$sigma
    limma_t <- limma_result$bayes$t %>%
        as.data.frame
    limma_p_adjust <- as_tibble(limma_pval) %>%
        mutate(across(everything(), ~p.adjust(.x, method = "BH")))

    limma_pval_tb <- setNames(limma_pval, paste0("pval", colnames(limma_pval)))
    limma_beta_tb <- setNames(limma_beta, paste0("beta", colnames(limma_beta)))
    limma_sigma_tb <- setNames(tibble(limma_sigma), c(str_glue("sigma{colnames(limma_result$fit$contrast)}")))
    limma_t_tb <- setNames(limma_t, paste0("tstat", colnames(limma_t)))
    limma_p_adjust_tb <- setNames(limma_p_adjust, paste0("fdr.pval", colnames(limma_pval)))

    limma_coef_tb <- bind_cols(limma_pval_tb, limma_beta_tb, limma_sigma_tb, limma_t_tb, limma_p_adjust_tb)
    limma_coef_tb <- limma_coef_tb %>% 
        mutate(rowname = rownames(limma_pval_tb)) %>%
        select(rowname, everything())
    return(limma_coef_tb)
}

limma_dx_meta_fit <- function(
        library_meta, library_class_counts,
        library_name = "library_id", count_name = "subcluster_ct", class_name = "ct_subcluster", 
        formula, randomize = FALSE, seed = 0, n = 1
    ) {

    library_class_mat <- pivot_matrix(library_class_counts, library_name, count_name, class_name, values_fill = 0)
    library_dx_other <- library_meta
    limma_design <- model.matrix(as.formula(formula), data = library_meta)
    limma_contrasts <- makeContrasts(
        contrasts = c(
            "sexM",
            "pmi",
            "age",
            "mean_percent_mito",
            "median_genes"
        ),
        levels = colnames(limma_design)
    )
    limma_result <- limma_voom_norm(library_class_mat, limma_design, limma_contrasts) 
    limma_result$form <- formula

    return(limma_result)
}


# Get permuted / randomzied p-values by comparing base t-statistic vs. permuted t-statistic
# for test dx{dx} - dxOther.
# if original t-stat is negative, the permuted p-value is the proportion of permuted t-values less than original.
# if original t-stat is positive, ... greater than original.
# additionally do p.adjust on permuted p-value within region / celltype.
summarize_dx <- function(base, permute, dx) {
    t_value_colname <- str_glue("tstatdx{dx} - dxOther")
    pval_tb <- base %>%
        group_by(region, cell_type, ct_subcluster) %>%
        group_nest(keep = TRUE) %>%
        mutate(pval_randomized_t = map_dbl(data, function(data) {
            permute_subset <- permute %>%
                filter(region %in% data$region, ct_subcluster %in% data$ct_subcluster)
            orig_t <- data[[t_value_colname]]
            if (orig_t < 0)
                n_gt <- sum(permute_subset[[t_value_colname]] < orig_t)
            else
                n_gt <- sum(permute_subset[[t_value_colname]] > orig_t)
            return(n_gt / nrow(permute_subset))
        })) %>%
        select(-data)
    pval_tb %>%
        group_by(region, cell_type) %>%
        mutate(fdr_pval_randomized_t = p.adjust(pval_randomized_t, method = "BH"))
}

write_limma <- function(tb, path) {
    tb %>%
        unnest(c(limma_coefs, form)) %>%
        select(-data, -limma_results) %>%
        write_csv(path)
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from, ...) {
    stopifnot(c(cols_from, values_from, rows_from) %in% colnames(tb))
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = all_of(values_from), ...) %>%
        ungroup
    
    tb_matrix <- select(tb_pivot, -all_of(rows_from)) %>% as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(all_of(rows_from)), collapse = "|"), .groups = "drop") %>%
        pluck("join_rowname")

    return(tb_matrix)
}

if (!interactive()) { main() }
