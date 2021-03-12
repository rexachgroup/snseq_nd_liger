# heatmap of top up / down lme per dx.
liblist <- c("tidyverse", "ComplexHeatmap", "png", "grid", "gridExtra")
l <- lapply(liblist, require, character.only = TRUE)

dge_list <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_wk.rds"
jacc_out <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_lme_jaccard.csv"
jacc_plot <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_lme_jaccard.pdf"

ngenes <- 150

main <- function() {
    subcluster_dges <- readRDS(dge_list)
    subcluster_filter <- subcluster_dges %>%
        filter(cluster_cell_type == "microglia")
    subcluster_cmp <- combn(as.character(subcluster_filter$ct_subcluster), 2) %>%
        t %>%
        as.data.frame
    colnames(subcluster_cmp) <- c("cmp_1", "cmp_2")
    subcluster_cmp <- as_tibble(subcluster_cmp)
    filter_args <- cross_df(
        list(dx_col = c("clinical_dxAD.estimate","clinical_dxPSP-S.estimate","clinical_dxbvFTD.estimate"), 
             type = c("up", "down"))
    )
    # full join to get comparison df
    cmp_df <- full_join(filter_args, subcluster_cmp, by = character())
    # do comparison on top n genes per subcluster
    cmp_df <- cmp_df %>%
        mutate(jacc_index = pmap_dbl(cmp_df, function(...) {
            current_row <- list(...)
            print(str_glue("{current_row$cmp_1} {current_row$cmp_2} {current_row$dx_col} {current_row$type}"))
            genes_1 <- subcluster_filter %>% 
                filter(ct_subcluster == current_row$cmp_1) %>% 
                pluck("broom_join", 1) %>%
                filter_dx_direction_estimate(., current_row$dx_col, current_row$type) %>%
                slice_head(n = ngenes) %>%
                pluck("gene")
            genes_2 <- subcluster_filter %>%
                filter(ct_subcluster == current_row$cmp_2) %>% 
                pluck("broom_join", 1) %>%
                filter_dx_direction_estimate(., current_row$dx_col, current_row$type) %>%
                slice_head(n = ngenes) %>%
                pluck("gene")

            return(jacc(genes_1, genes_2))
        }))

    write_csv(cmp_df, jacc_out) 

    # Iterate over comparisons and create jaccard plot for each.
    jacc_heatmaps <- pmap(filter_args, function(...) {
        cr <- list(...)
        print(str_glue("{ngenes} {cr$type} {cr$dx_col}"))
        cmp_mat <- cmp_df %>%
            filter(dx_col == cr$dx_col, type == cr$type) %>%
            select(cmp_1, cmp_2, jacc_index) %>%
            pivot_wider(names_from = "cmp_2", values_from = "jacc_index") %>%
            as.data.frame()
        rownames(cmp_mat) <- cmp_mat[, "cmp_1"]
        cmp_mat <- cmp_mat[, -1]
        cmp_mat <- data.matrix(cmp_mat)
        print(str_glue("{dim(cmp_mat)}"))

        jacc_title <- str_glue("jaccard overlap for top {cr$type} {ngenes} ranked by {cr$dx_col}")
        jacc_heatmap <- Heatmap(cmp_mat, cluster_rows = FALSE, cluster_columns = FALSE, name = "jaccard", column_title = jacc_title)
        jacc_rastr <- tempfile()
        print(jacc_title)
        png(jacc_rastr, width = 10, height = 10, units = "in", res = 200)
        print(jacc_heatmap)
        graphics.off()

        rasterGrob(readPNG(jacc_rastr))
    })

    ggsave(jacc_plot, marrangeGrob(jacc_heatmaps, newpage = FALSE, nrow = 1, ncol = 1), width = 10, height = 10)
}


filter_dx_direction_estimate <- function(tb, col_selector = NULL, type = NULL) {
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
