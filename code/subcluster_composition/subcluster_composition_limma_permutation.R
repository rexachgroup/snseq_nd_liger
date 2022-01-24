# subcluster composition permutation testing.
liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl", "future.apply", "progressr")
lapply(liblist, require, character.only = TRUE)
set.seed(NULL)
handlers("txtprogressbar")

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma_permutation/"
AZIMUTH_CLUSTERS <- "../../ext/azimuth/data/meta.csv"

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
        filter(cell_type == "excitatory") %>%
        group_by(region, library_id, ct_subcluster) %>%
        summarize(
            cell_type = unique(cell_type), 
            liger_cluster = unique(liger_clusters),
            subcluster_ct = n()
        ) %>%
        group_by(region, cell_type) %>%
        group_nest(keep = TRUE)

    # Get nominal model fit.
    library_subcluster_other_limma <- library_subcluster_counts %>%
        mutate(limma_other_results = map(data, function(.x) {
            library_subset_meta <- filter(library_meta, library_id %in% .x$library_id)
            fit <- limma_dx_other_fit(
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

    # Permute 10000x with randomized labels.
    plan(multicore)
    library_subcluster_permute_limma <- library_subcluster_counts %>%
        mutate(limma_other_results = map(data, function(.x) {
            library_subset_meta <- filter(library_meta, library_id %in% .x$library_id)
            fit <- limma_dx_other_fit(
                library_subset_meta,
                .x,
                library_name = "library_id",
                count_name = "subcluster_ct",
                class_name = "ct_subcluster",
                formula = "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes",
                randomize = TRUE,
                n = 10000
            )
            return(fit)
        }))

    saveRDS(library_subcluster_permute_limma, file.path(OUT_DIR, "subcluster_composition_randomized.rds"))

    # Extract limma result per iteration and join into large table.
    library_subcluster_permute_limma <- library_subcluster_permute_limma %>%
        mutate(limma_result_fmt = map(limma_other_results, function(.x) {
            cross_dx_list <- pmap(.x, function(...) {
                cr <- list(...)
                writeLines(str_glue("{cr$test_dx}"))
                progressr::with_progress({
                    progressor <- progressr::progressor(steps = length(cr$limma_result))
                    limma_coefs <- future_Map(function(i) {
                        progressor(i)
                        coefs <- as_tibble(limma_extract_coefs(cr$limma_result[[i]]))
                        rownames(coefs) <- NULL
                        return(coefs)
                    }, seq_along(cr$limma_result))
                })
                return(bind_rows(limma_coefs))
            })
            cross_dx_tb <- cross_dx_list %>%
                bind_cols
            rowname <- cross_dx_tb[["rowname...1"]]
            cross_dx_tb <- cross_dx_tb %>%
                select(-contains("rowname")) %>%
                mutate(rowname = rowname) %>%
                select(rowname, everything())
            return(cross_dx_tb)
        }))

    library_subcluster_permute_tb <- library_subcluster_permute_limma %>%
        select(-data, -limma_other_results) %>%
        unnest(limma_result_fmt) %>%
        rename(ct_subcluster = rowname)

    # same for original run -- extract coefficients per dx test and join into large table.
    library_subcluster_other_tb <- library_subcluster_other_limma %>%
        mutate(limma_result_fmt = map(limma_other_results, function(limma_result_tb) {
            cross_dx <- pmap(limma_result_tb, function(...) {
                cr <- list(...)
                writeLines(str_glue("{cr$test_dx}"))
                dx_tb <- limma_extract_coefs(cr$limma_result[[1]])
            })
            cross_dx_tb <- cross_dx %>%
                map(~select(., -rowname)) %>%
                bind_cols %>%
                rownames_to_column
            return(cross_dx_tb)
        })) %>%
        unnest(limma_result_fmt) %>%
        select(-data, -limma_other_results) %>%
        rename(ct_subcluster = rowname)
    
    write_csv(library_subcluster_other_tb, file.path(OUT_DIR, "subcluster_composition_dx_other.csv"))
    write_csv(library_subcluster_permute_tb, file.path(OUT_DIR, "subcluster_composition_dx_permute.csv"))

    # get permutation test p-value per dx.
    subcluster_ad <- summarize_dx(library_subcluster_other_tb, library_subcluster_permute_tb, "AD")
    subcluster_psp <- summarize_dx(library_subcluster_other_tb, library_subcluster_permute_tb, "PSP.S")
    subcluster_ftd <- summarize_dx(library_subcluster_other_tb, library_subcluster_permute_tb, "bvFTD")
    write_csv(subcluster_ad, file.path(OUT_DIR, "subcluster_composition_dx_permute_ad.csv"))
    write_csv(subcluster_psp, file.path(OUT_DIR, "subcluster_composition_dx_permute_psp.csv"))
    write_csv(subcluster_ftd, file.path(OUT_DIR, "subcluster_composition_dx_permute_ftd.csv"))

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

limma_dx_other_fit <- function(
        library_meta, library_class_counts,
        library_name = "library_id", count_name = "subcluster_ct", class_name = "ct_subcluster", 
        formula, randomize = FALSE, seed = 0, n = 1
    ) {

    disease_dx <- c("AD", "bvFTD", "PSP.S")
    library_class_mat <- pivot_matrix(library_class_counts, library_name, count_name, class_name, values_fill = 0)
    # cross join with disease dx to get 3x meta tables
    library_dx_other <- library_meta %>% 
        full_join(tibble(test_dx = disease_dx), by = character()) %>%
        group_by(test_dx) %>%
        group_nest(keep = TRUE)

    prev_seed <- .Random.seed
    set.seed(seed)
    limma_dx_other <- library_dx_other %>%
        mutate(limma_result = pmap(., function(...) {
            cr <- list(...)
            future_Map(function(i) {
                data <- cr$data %>%
                    mutate(dx = dx_other_relabel(dx, test_lvl = cr$test_dx, randomize = randomize))

                if (all(data$dx == "Other")) {
                    return(NA)
                }
                limma_design <- model.matrix(as.formula(formula), data = data)
                contrast_str <- as.character(str_glue("dx{unique(data$test_dx)} - dxOther"))
                limma_contrasts <- makeContrasts(
                    contrasts = c(contrast_str),
                    levels = colnames(limma_design)
                )
                limma_result <- limma_voom_norm(library_class_mat, limma_design, limma_contrasts) 
                limma_result$form <- formula
                return(limma_result)
            }, 1:n, future.seed = seed)
        }))
    .Random.seed <- prev_seed
    return(limma_dx_other)
}


# relabel dx, with background levels set to Other.
# if randomize = true, shuffle labels.
dx_other_relabel <- function(lvl_vec, test_lvl, randomize = TRUE) {
    if(is.character(lvl_vec))
        lvl_vec <- as.factor(lvl_vec)
    bg_levels  <- setdiff(levels(lvl_vec), test_lvl)
    if (randomize) {
        order_rand <- sample(lvl_vec, size = length(lvl_vec), replace = FALSE)
    } else {
        order_rand <- lvl_vec
    }
    rand_dx_vec <- fct_other(order_rand, drop = bg_levels)
    return(rand_dx_vec)
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
