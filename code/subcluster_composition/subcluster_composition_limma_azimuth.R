# subcluster composition calculation and run on ext/Azimuth excitatory clusters.
liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma_azimuth/"
AZIMUTH_CLUSTERS <- "../../ext/azimuth/data/meta.csv"

main <- function() {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    azimuth_meta <- read_csv(AZIMUTH_CLUSTERS)

    meta <- azimuth_meta %>%
        mutate(
            clinical_dx = fct_relevel(clinical_dx, "Control"),
            predicted_cluster_celltype = str_extract(predicted.cluster, "^\\w+")
        )

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

    # filter by predicted clusters that have control samples.
    ctl_filter <- meta %>% group_by(predicted.cluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)

    # Count the number of cells belonging to each sample for each predicted cluster.
    library_cluster_counts <- meta %>%
        filter(predicted.cluster %in% ctl_filter$predicted.cluster,
               predicted.cluster.score > 0.5) %>%
        group_by(region, library_id, predicted.cluster) %>%
        filter(predicted_cluster_celltype == "Exc") %>%
        summarize(
            cluster_ct = n(),
            predicted_cluster_celltype = unique(predicted_cluster_celltype),
            .groups = "drop"
        )  %>%
        group_by(region, predicted_cluster_celltype) %>%
        group_nest(keep = TRUE)

    library_cluster_limma <- library_cluster_counts %>%
        mutate(limma_results = map(data, function(.x) {
            library_subset_meta <- filter(library_meta, library_id %in% .x$library_id)
            tryCatch(limma_dx_contrast_fit(
                library_subset_meta,
                .x,
                "predicted.cluster",
                "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes"
            ), error = function(x) return(NA))
        }))
    limma_fit <- library_cluster_limma %>%
        mutate(
            limma_coefs = map(limma_results, function(.x) {
                coef_tb <- limma_extract_coefs(.x)
            }),
            form = map(limma_results, function(.x) { .x$form })
        )


    write_limma(limma_fit, file.path(OUT_DIR, "cluster_composition_azimuth_exc.tsv"))

    
    library_subclass_counts <- meta %>%
        filter(predicted.subclass.score > 0.5) %>%
        group_by(region, library_id, predicted.subclass) %>%
        summarize(cluster_ct = n()) %>%
        group_by(region) %>%
        group_nest(keep = TRUE)

    library_subclass_limma <- library_subclass_counts %>%
        mutate(limma_results = map(data, function(.x) {
            library_subset_meta <- filter(library_meta, library_id %in% .x$library_id)
            return(limma_dx_contrast_fit(
                library_subset_meta,
                .x,
                "predicted.subclass",
                "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes"
            ))
        }))

    library_subclass_limma <- library_subclass_limma %>%
        mutate(
            limma_coefs = map(limma_results, function(.x) {
                coef_tb <- limma_extract_coefs(.x)
            }),
            form = map(limma_results, function(.x) { .x$form })
        )

    write_limma(library_subclass_limma, file.path(OUT_DIR, "subclass_composition_azimuth.tsv"))
    

}

limma_dx_contrast_fit <- function(library_meta, library_class_counts, class_name, formula) {
    library_class_mat <- pivot_matrix(library_class_counts, "library_id", "cluster_ct", class_name, values_fill = 0)
    limma_design <- model.matrix(as.formula(formula), data = library_meta)
    limma_contrasts <- makeContrasts(
        contrasts = c(
            "dxAD - dxControl",
            "dxbvFTD - dxControl",
            "dxPSP.S - dxControl",
            "(dxAD + dxbvFTD + dxPSP.S)/3 - dxControl",
            "dxAD - (dxControl + dxbvFTD + dxPSP.S) / 3",
            "dxbvFTD - (dxControl + dxAD + dxPSP.S) / 3",
            "dxPSP.S - (dxControl + dxAD + dxbvFTD) / 3",
            "dxbvFTD - (dxAD + dxPSP.S) / 2"
        ),
        levels = colnames(limma_design)
    )
    limma_result <- limma_voom_norm(library_class_mat, limma_design, limma_contrasts)
    limma_result$form <- formula
    return(limma_result)
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
    limma_p_adjust <- as_tibble(limma_pval) %>%
        mutate(across(everything(), ~p.adjust(.x, method = "BH")))

    limma_pval_tb <- setNames(limma_pval, paste0("pval", colnames(limma_pval)))
    limma_beta_tb <- setNames(limma_beta, paste0("beta", colnames(limma_beta)))
    limma_sigma_tb <- setNames(tibble(limma_sigma), c(str_glue("sigma{colnames(limma_result$fit$contrast)}")))
    limma_p_adjust_tb <- setNames(limma_p_adjust, paste0("fdr.pval", colnames(limma_pval)))

    limma_coef_tb <- bind_cols(limma_pval_tb, limma_beta_tb, limma_sigma_tb, limma_p_adjust_tb)
    limma_coef_tb <- limma_coef_tb %>% 
        mutate(rowname = rownames(limma_pval_tb)) %>%
        select(rowname, everything())
    return(limma_coef_tb)
}

write_limma <- function(tb, path) {
    tb %>%
        unnest(c(limma_coefs, form)) %>%
        select(-data, -limma_results) %>%
        write_csv(path)
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from, ...) {
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
