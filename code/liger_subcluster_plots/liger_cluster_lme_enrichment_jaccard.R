# heatmap of top up/down for gene signatures lme
liblist <- c("tidyverse", "ComplexHeatmap", "png", "grid", "gridExtra", "readxl")
l <- lapply(liblist, require, character.only = TRUE)

dge_list <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
jacc_out <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_enriched_lme_jaccard.csv"
jacc_plot <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_enriched_lme_jaccard.pdf"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"

main <- function() {
    filter_args <- tribble(
        ~cluster_cell_type, ~fdr_max_cutoff, ~estimate_min_cutoff,
        "excitatory", 0.05, 0.2,
        "inhibitory", 0.05, 0.2,
        "oligodendrocyte", 0.05, 0.2,
        "astrocyte", 0.05, 0.2,
        "microglia", 0.1, 0.1,
        "opc", 0.1, 0.1,
        "endothelia", 0.1, 0.1
    )
    excludes <- read_xlsx(clusters_exclude_file)
    subcluster_dges <- readRDS(dge_list) %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster)

    cmp_df <- pmap(filter_args, function(...) {
        cr_f <- list(...)
        print(str_glue("{cr_f$cluster_cell_type} {cr_f$fdr_max_cutoff} {cr_f$estimate_min_cutoff}"))
        subcluster_filter <- subcluster_dges %>%
            filter(cluster_cell_type %in% cr_f$cluster_cell_type)
        subcluster_cmp <- cross_df(list(cmp_1 = subcluster_filter$ct_subcluster, cmp_2 = subcluster_filter$ct_subcluster))
        subcluster_cmp %>%
            mutate(jacc_index = pmap_dbl(subcluster_cmp, function(...) {
                cr <- list(...)
                #print(str_glue("{cr$cmp_1} {cr$cmp_2}"))
                genes_1 <- subcluster_filter %>%
                    filter(ct_subcluster == cr$cmp_1) %>% 
                    pluck("broom_join", 1) %>%
                    filter(across(contains("p.value.adj")) < cr_f$fdr_max_cutoff, abs(across(contains("estimate"))) > cr_f$estimate_min_cutoff) %>%
                    pluck("gene")
                genes_2 <- subcluster_filter %>%
                    filter(ct_subcluster == cr$cmp_2) %>% 
                    pluck("broom_join", 1) %>%
                    filter(across(contains("p.value.adj")) < cr_f$fdr_max_cutoff, abs(across(contains("estimate"))) > cr_f$estimate_min_cutoff) %>%
                    pluck("gene")
                return(jacc(genes_1, genes_2))
            }))
    }) %>%
    setNames(nm = filter_args$cluster_cell_type)

    write_csv(bind_rows(cmp_df), jacc_out) 

    # Iterate over comparisons and create jaccard plot for each.
    jacc_heatmaps <- pmap(filter_args, function(...) {
        cr <- list(...)
        cmp_mat <- cmp_df[[cr$cluster_cell_type]] %>%
            select(cmp_1, cmp_2, jacc_index) %>%
            pivot_wider(names_from = "cmp_2", values_from = "jacc_index") %>%
            as.data.frame()
        rownames(cmp_mat) <- cmp_mat[, "cmp_1"]
        cmp_mat <- cmp_mat[, -1]
        cmp_mat <- data.matrix(cmp_mat)
        print(str_glue("{dim(cmp_mat)}"))

        jacc_title <- str_glue(
            "jaccard overlap for {cr$cluster_cell_type} clusters", 
            "(fdr < {cr$fdr_max_cutoff}, abs(estimate) > {cr$estimate_min_cutoff})",
            "per liger subcluster vs. all-celltype background", .sep = "\n")
        jacc_heatmap <- Heatmap(cmp_mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "jaccard", column_title = jacc_title)
        #jacc_rastr <- tempfile()
        print(jacc_title)
        #png(jacc_rastr, width = 10, height = 10, units = "in", res = 200)
        #print(jacc_heatmap)
        #graphics.off()

        #rasterGrob(readPNG(jacc_rastr))
        return(jacc_heatmap)
    })

    pdf(jacc_plot, width = 10, height = 10)
    walk(jacc_heatmaps, draw)
    dev.off()
}


filter_direction_estimate <- function(tb, col_selector = NULL, type = NULL) {
    valid_types = c("up", "down")
    stopifnot(type %in% valid_types)

    if (type == "up") {
        return(tb %>% 
               arrange(desc(across(contains(col_selector)))) %>% 
               filter(across(contains(col_selector)) > 0))
    } else if (type == "down") { 
        return(tb %>% 
               arrange(across(contains(col_selector))) %>% 
               filter(across(contains(col_selector)) < 0))
    }
}

jacc <- function(x, y) {
    return(length(intersect(x, y)) / length(union(x, y)))
}

if (!interactive()) {
    main()
}
