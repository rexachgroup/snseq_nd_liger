# limma across dx contrasts of cluster counts with haplotype.
liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/haplotype/"
SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
SUBJ_META <- "../../resources/20220119 SampleLevel_metadata_2021_rexach_v3_brainbank_LL_PIDN ApoE and Tau_aln20220225.xlsx"

main <- function() {
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    liger_meta <- readRDS(META_FILE)
    subj_meta <- read_xlsx(SUBJ_META) %>%
        filter(!is.na(Tau_H)) %>%
        mutate(autopsy_id = paste0("P", PIDN), library_id = library_ID...6) %>%
        select(library_id, Tau_H)
    
    liger_meta_hap <- liger_meta %>%
        inner_join(subj_meta, by = "library_id") %>%
        filter(Tau_H %in% c("H1/H1", "H1/H2"))
    
    ctl_filter <- liger_meta %>% group_by(ct_subcluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)

    library_meta <- liger_meta_hap %>%
        group_by(library_id) %>%
        summarize(
            dx = make.names(unique(clinical_dx)),
            region = unique(region),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            tau = make.names(unique(Tau_H))
        )
    
    library_subcluster_counts <- liger_meta_hap %>%
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
                formula = "~0 + tau + dx + pmi + age + sex + mean_percent_mito + median_genes",
                randomize = FALSE
            )
            return(fit)
        }))
    
    library_subcluster_meta_tb <- library_subcluster_meta_limma %>%
        mutate(coefs = map(limma_meta_results, limma_extract_coefs)) %>%
        select(-data, -limma_meta_results) %>%
        unnest(coefs)
    write_csv(library_subcluster_meta_tb, file.path(OUT_DIR, "subcluster_composition_haplotype.csv"))
}

limma_voom_norm <- function(data, design, contrast) {
    dgelist <- DGEList(data)
    dgelist_norm <- calcNormFactors(dgelist, method = "TMM")
    elist_voom <- voom(dgelist_norm, design)
    dx_fit <- lmFit(elist_voom, design)
    dx_contrast_fit <- contrasts.fit(dx_fit, contrast)
    bayes <- tryCatch(eBayes(dx_contrast_fit), error = function(x) NA)
    return(list(fit = dx_contrast_fit, bayes = bayes))
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

limma_dx_meta_fit <- function(
        library_meta, library_class_counts,
        library_name = "library_id", count_name = "subcluster_ct", class_name = "ct_subcluster", 
        formula, randomize = FALSE, seed = 0, n = 1
    ) {

    library_class_mat <- pivot_matrix(library_class_counts, library_name, count_name, class_name, values_fill = 0)
    limma_design <- model.matrix(as.formula(formula), data = library_meta)
    limma_contrasts <- makeContrasts(
        contrasts = c(
            "tauH1.H2 - tauH1.H1",
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
