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
              subcluster_ct = n())

library_celltype_counts_full <- meta %>%
    filter(cluster_cell_type == cell_type) %>%
    group_by(region, library_id, cell_type) %>%
    summarize(
        clinical_dx = unique(clinical_dx),
        age = unique(age),
        pmi = unique(pmi),
        pmi = unique(pmi),
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
            dx_bayes <- eBayes(limma)
            dx_bayes$p.value %>%
                as_tibble(rownames = "ct_subcluster")
        }),
        formula = form
    )
}

summarize_limma <- function(tb) {
    tb %>%
        unnest(limma_pval) %>%
        group_by(region) %>%
        summarize(across(matches("dx"), function(x) sum(!is.na(x) & x < 0.05))) %>%
        print(width = Inf)
}

write_limma <- function(tb, path) {
    tb %>% 
        unnest(limma_pval) %>%
        select(-data, -matrix, -limma) %>%
        write_tsv(path)
}


limma_dx <- limma_mutate(library_cluster_mats, form = "~0 + dx")
limma_dx_log_pmi <- limma_mutate(library_cluster_mats, form = "~0 + dx + log_pmi + age + mean_percent_mito + median_genes")
limma_dx_pmi <- limma_mutate(library_cluster_mats, form = "~0 + dx + pmi + age + mean_percent_mito + median_genes")

summarize_limma(limma_dx)
summarize_limma(limma_dx_log_pmi)
summarize_limma(limma_dx_pmi)

write_limma(limma_dx, file.path(OUT_DIR, "subcluster_composition_dx.tsv"))
write_limma(limma_dx_log_pmi, file.path(OUT_DIR, "subcluster_composition_dx_log_pmi.tsv"))
write_limma(limma_dx_pmi, file.path(OUT_DIR, "subcluster_composition_dx_pmi.tsv"))

saveRDS(limma_dx, file.path(OUT_DIR, "subcluster_composition_dx.rds"))
saveRDS(limma_dx_log_pmi, file.path(OUT_DIR, "subcluster_composition_dx_log_pmi.rds"))
saveRDS(limma_dx_pmi, file.path(OUT_DIR, "subcluster_composition_dx_pmi.rds"))

