# limma across dx contrasts / dx splits within a celltype.
liblist <- c("tidyverse", "broom", "edgeR", "limma")
lapply(liblist, require, character.only = TRUE)

META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/limma/"
SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

meta_tb <- readRDS(META_FILE)
percluster_tb <- readRDS(SUBCLUSTER_FILTER_FILE) %>% 
    filter(subcluster_3libs_above_10umi) %>%
    mutate(ct_subcluster = paste(region, ct_subcluster, sep = "-")) %>%
    mutate(ct_subcluster = fct_collapse(ct_subcluster, "calcarine-excitatory-210" = c("calcarine-excitatory-2", "calcarine-excitatory-10")))

# lump together calcarine-excitatory-2 and calcarine-excitatory-10
meta <- meta_tb %>%
    mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
        clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
    mutate(ct_subcluster = fct_collapse(ct_subcluster, "calcarine-excitatory-210" = c("calcarine-excitatory-2", "calcarine-excitatory-10")),
        liger_clusters = ifelse((region == "calcarine" & cluster_cell_type == "excitatory"), 210, liger_clusters)
    )

ctl_filter <- meta %>% group_by(ct_subcluster) %>% summarize(ctl = any(clinical_dx == "Control")) %>% filter(ctl == TRUE)

# Count the number of cells belonging to a control sample for each subcluster.
library_subcluster_counts <- meta %>%
    filter(cluster_cell_type == cell_type,
           ct_subcluster %in% percluster_tb$ct_subcluster,
           ct_subcluster %in% ctl_filter$ct_subcluster,) %>%
    group_by(region, library_id, ct_subcluster) %>%
    summarize(cell_type = unique(cell_type), 
              liger_cluster = unique(liger_clusters),
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



# for each data + matrix row, fit limma model using form, then compute contrasts and estimate eBayes
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
                "dxPSP.S - (dxControl + dxAD + dxbvFTD) / 3"
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
    tb %>%
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
        })
    )
}

summarize_limma <- function(tb) {
    tb %>%
        filter(!is.na(limma_pval)) %>%
        mutate(limma_cmb = pmap(list(limma_beta, limma_pval, limma_sigma), 
            function(beta, pval, s2) {
                inner_join(beta, pval, by = "ct_subcluster") %>%
                    inner_join(s2, by = "ct_subcluster")
            })
        ) %>%
        unnest(limma_cmb)
}

write_limma <- function(tb, path) {
    tb %>% 
        summarize_limma %>%
        select(-data, -matrix, -limma, -limma_bayes, -limma_pval, -limma_beta, -limma_sigma) %>%
        write_tsv(path)
}

limma_dx_pmi <- limma_fit_contrast(library_cluster_mats, 
    form = "~0 + dx + pmi + age + sex + mean_percent_mito + median_genes")
limma_dx_pmi_coefs <- limma_extract_coefs(limma_dx_pmi)
summarize_limma(limma_dx_pmi_coefs)

write_limma(limma_dx_pmi_coefs, file.path(OUT_DIR, "subcluster_composition_dx_pmi.tsv"))
saveRDS(limma_dx_pmi_coefs, file.path(OUT_DIR, "subcluster_composition_dx_pmi.rds"))

