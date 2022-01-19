# limma across dx contrasts / dx splits within specific subclusters.
liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma_subclusters/"
SUBCLUSTER_FILTER_FILE <- "../../resources/subclusters_removed_byQC_final.xlsx"

main <- function() {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    meta_tb <- readRDS(META_FILE)
    excludes <- read_xlsx(SUBCLUSTER_FILTER_FILE)

    meta <- meta_tb %>%
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
        mutate(ct_subcluster = fct_collapse(ct_subcluster, "insula-inhibitory-0/6" = c("insula-inhibitory-0", "insula-inhibitory-6")))

    ctl_filter <- meta %>% group_by(ct_subcluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)

    library_subcluster_counts <- meta %>%
        filter(cluster_cell_type == cell_type,
               !ct_subcluster %in% excludes$ct_subcluster,
               ct_subcluster %in% ctl_filter$ct_subcluster) %>%
        filter(region == "insula", cell_type == "inhibitory") %>%
        group_by(region, library_id, ct_subcluster) %>%
        summarize(cell_type = unique(cell_type), 
                  subcluster_ct = n())

    # Count the total number of cells per Seurat cluster.
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

    library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type"))
    
    library_subset_mats <- library_pct %>%
        group_by(region, cell_type) %>%
        group_nest(keep = TRUE) %>%
        mutate(matrix = map(data, function(data) {
            pivot_matrix(data, "library_id", "subcluster_ct", "ct_subcluster", values_fn = sum, values_fill = 0)
        }))
    
    form <- "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes"
    mes_fit <- limma_fit_contrast(library_subset_mats, form = form) %>%
        limma_extract_coefs

    write_limma(mes_fit, file.path(OUT_DIR, "subcluster_composition_mes.tsv"))
    saveRDS(mes_fit, file.path(OUT_DIR, "subcluster_composition_mes.rds"))
}

# for each data + matrix row, fit limma model using form, then compute contrasts and estimate eBayes
limma_fit_contrast <- function(tb, form) {
    tb %>%
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
            design <- model.matrix(as.formula(form), data = dx_fct)
            
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
limma_extract_coefs <- function(tb) {
    tb_fmt <- tb %>%
    mutate(
        limma_pval = map(limma_bayes, function(dx_bayes) {
            dx_bayes$p.value %>%
                as.data.frame %>%
                rename_with(function(x) paste0("pval", x), everything()) %>%
                as_tibble(rownames = "ct_subcluster")
        }),
        limma_beta = map(limma, function(limma) {
            limma$coefficients %>%
                as.data.frame %>%
                rename_with(function(x) paste0("beta", x), everything()) %>%
                as_tibble(rownames = "ct_subcluster")
        }),
        limma_sigma = map(limma_bayes, function(dx_bayes) {
            tibble(ct_subcluster = rownames(dx_bayes$p.value), sigma = dx_bayes$sigma)
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
            Reduce(function(.x, .y) inner_join(.x, .y, by = "ct_subcluster"), cr)
        })) %>%
        unnest(limma_cmb)
}

write_limma <- function(tb, path) {
    tb %>% 
        summarize_limma %>%
        select(-data, -matrix, -limma, -limma_bayes, -limma_pval, -limma_beta, -limma_sigma, -limma_p_adjust) %>%
        write_tsv(path)
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
