# Copy subcluster_composition_limma. Merge insula-excitatory 2+6 and precg-excitatory 4+5.

liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl", "ggpubr")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma_merge/"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"

main <- function() {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    meta_tb <- readRDS(META_FILE)
    excludes <- read_xlsx(clusters_exclude_file)
    percluster_tb <- readRDS(SUBCLUSTER_FILTER_FILE) %>%
        filter(subcluster_3libs_above_10umi) %>%
        mutate(ct_subcluster = paste(region, ct_subcluster, sep = "-"))

    meta <- meta_tb %>%
        filter(cluster_cell_type == "excitatory") %>%
        mutate(
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            clinical_dx = fct_relevel(clinical_dx, "Control")
        ) %>%
        mutate(
            ct_subcluster = fct_drop(fct_relevel(fct_recode(ct_subcluster,
                "insula-excitatory-2-6" = "insula-excitatory-2",
                "insula-excitatory-2-6" = "insula-excitatory-6",
                "precg-excitatory-4-5" = "preCG-excitatory-4",
                "precg-excitatory-4-5" = "preCG-excitatory-6",
            ), ~str_sort(., numeric = T)))
        )
        

    ctl_filter <- meta %>% 
        group_by(ct_subcluster) %>% 
        summarize(ctl = any(clinical_dx == "Control")) %>% 
        filter(ctl == TRUE)

    # Count the number of cells for each subcluster that has control samples.
    library_subcluster_counts <- meta %>%
        filter(cluster_cell_type == cell_type,
               ct_subcluster %in% ctl_filter$ct_subcluster) %>%
        group_by(region, library_id, ct_subcluster) %>%
        summarize(
            cell_type = unique(cell_type), 
            subcluster_ct = n(),
            .groups = "keep"
        )

    # Count the total number of cells + summarize metadata per library + Seurat cluster subset.
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
            .groups = "keep"
        )

    # Join tables and pivot to library x subcluster cell counts matrix per Seurat cluster.
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

    limma_dx_pmi <- limma_fit_contrast(library_cluster_mats, 
        form = "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes")
    limma_dx_pmi_coefs <- limma_extract_coefs(limma_dx_pmi, excludes)
    summarize_limma(limma_dx_pmi_coefs)

    library_ct_patchwork <- library_pct %>%
        group_by(ct_subcluster) %>%
        group_split() %>%
        map(function(data) {
            ggplot_subcluster_library_counts(data) +
            ggtitle(str_glue("{unique(data$ct_subcluster)} subcluster count, control + umi filter"))
        })

    write_csv(library_pct, file.path(OUT_DIR, "subcluster_composition_ct.csv"))
    write_limma(limma_dx_pmi_coefs, file.path(OUT_DIR, "subcluster_composition_dx_pmi.tsv"))
    saveRDS(limma_dx_pmi_coefs, file.path(OUT_DIR, "subcluster_composition_dx_pmi.rds"))
    
    pdf(file.path(OUT_DIR, "counts_boxplot.pdf"), width = 10, height = 10)
    print(library_ct_patchwork)
    graphics.off()
}

# for each metadata + matrix row (one Seurat cluster), fit limma model using form, then compute contrasts and estimate eBayes
limma_fit_contrast <- function(tb, form) {
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
                "dxPSP.S - (dxControl + dxAD + dxbvFTD) / 3",
                "dxbvFTD - (dxAD + dxPSP.S) / 2"
                ),
                levels = colnames(coef(dx_fit))
            )
            dx_contrast_fit <- contrasts.fit(dx_fit, limma_contrasts)
            return(dx_contrast_fit)
        }),
        limma_bayes = map(limma, function(limma) { 
            tryCatch({
                dx_bayes <- eBayes(limma)
            }, error = function(x) NA)
        }),
        formula = form
    )
}

# extract + format eBayes coefficients
limma_extract_coefs <- function(tb, excludes) {
    tb_fmt <- tb %>%
    mutate(
        limma_pval = map(limma_bayes, function(dx_bayes) {
            dx_bayes$p.value %>%
                as.data.frame %>%
                rename_with(function(x) paste0("pval", x), everything()) %>%
                as_tibble(rownames = "ct_subcluster") %>%
                filter(!ct_subcluster %in% excludes$ct_subcluster)
        }),
        limma_beta = map(limma, function(limma) {
            limma$coefficients %>%
                as.data.frame %>%
                rename_with(function(x) paste0("beta", x), everything()) %>%
                as_tibble(rownames = "ct_subcluster") %>%
                filter(!ct_subcluster %in% excludes$ct_subcluster)
        }),
        limma_sigma = pmap(list(limma_pval, limma_bayes), function(pval_tb, dx_bayes) {
            ct_i <- match(pval_tb$ct_subcluster, rownames(dx_bayes$p.value))
            tibble(ct_subcluster = rownames(dx_bayes$p.value)[ct_i], sigma = dx_bayes$sigma[ct_i])
        }),
        limma_p_adjust = map(limma_pval, function(pval_tb) {
            mutate(pval_tb, across(
                starts_with("pval"), 
                ~p.adjust(.x, method = "BH"), 
                .names = "fdr.{.col}")
            ) %>%
            select(ct_subcluster,starts_with("fdr"))
        })
    )
}

summarize_limma <- function(tb) {
    tb %>%
        filter(!is.na(limma_pval)) %>%
        mutate(limma_cmb = pmap(list(limma_beta, limma_pval, limma_p_adjust, limma_sigma), function(...) {
            cr <- list(...)
            Reduce(function(.x, .y) full_join(.x, .y, by = "ct_subcluster"), cr)
        })) %>%
        unnest(limma_cmb)
}

write_limma <- function(tb, path) {
    tb %>% 
        summarize_limma %>%
        select(-data, -matrix, -limma, -limma_bayes, -limma_pval, -limma_beta, -limma_sigma, -limma_p_adjust) %>%
        write_tsv(path)
}

ggplot_subcluster_library_counts <- function(data, ct_col = "subcluster_ct") {
    comparisons <- list(c("Control", "AD"), c("Control", "bvFTD"), c("Control", "PSP-S"))
    max_dat <- quantile(data[[ct_col]], 0.95)
    range <- max(data[ct_col]) - min(data[[ct_col]])
    comparison_y <- c(max_dat, max_dat + (0.05 * range), max_dat + (0.10 * range))
    ggplot(data, aes_string(x = "clinical_dx", y = ct_col, color = "clinical_dx")) +
        geom_boxplot() +
        geom_jitter() +
        stat_compare_means(data = data,comparisons = comparisons, label.y = comparison_y)
}

main()
