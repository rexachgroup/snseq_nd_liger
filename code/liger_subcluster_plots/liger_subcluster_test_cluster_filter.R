# Comparison of changes in percent_mito and number_umi before and after filtering on cluster_cell_type == cell_type.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_seurat_rds <-
  "../../analysis/pci_import/pci_seurat.rds"
in_seurat_liger <-
  "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"

## Outputs
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"
dir.create(out_path_base, recursive = TRUE)

main <- function() {
    liger_meta <- readRDS(in_seurat_liger) %>%
        mutate( 
            liger_clusters = fct_inseq(liger_clusters),
            ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            log_number_umi = log(number_umi)
        ) %>%
        mutate(
            ct_subcluster = factor(ct_subcluster, levels = str_sort(unique(ct_subcluster), numeric = TRUE))
        )
 
    meta_subset <- liger_meta %>% 
        group_by(ct_subcluster)
    
    meta_subset_celltype_match <- liger_meta %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(ct_subcluster)
    
    meta_sum <- meta_subset %>%
        group_by(ct_subcluster) %>%
        summarize(n_before = n())

    meta_sum_celltype_match <- meta_subset_celltype_match %>%
        group_by(ct_subcluster) %>%
        summarize(n_after = n())

    meta_celltype_counts <- map(unique(meta_subset$cell_type), function(ct) {
            meta_subset %>%
                group_by(ct_subcluster) %>%
                summarize(across(.cols = cell_type, .fns = ~length(which(.x == ct)), .names = "n_{ct}"))
        }) %>%
        reduce(inner_join, by = "ct_subcluster")
    
    meta_sum_join <- reduce(list(meta_sum, meta_sum_celltype_match, meta_celltype_counts), inner_join, by = "ct_subcluster") %>%
        mutate(across(.cols = -c(ct_subcluster, n_before), .fns = ~ .x / n_before, .names = "ratio_{.col}"))
    write_csv(meta_sum_join, file.path(out_path_base, "meta_sum.csv"))
    
    meta_merge <- inner_join(
        group_nest(meta_subset, keep = TRUE),
        group_nest(meta_subset_celltype_match, keep = TRUE),
        by = "ct_subcluster")


    comp_vars <- c("percent_mito", "number_umi")
    comp_plots <- pmap(meta_merge, function(...) {
            cur_row = list(...)
            var_plots <- map(comp_vars, function(comp_var) {
                aes_plot <- aes_string(x = "ct_subcluster", y = comp_var)
                unf <- ggplot(cur_row$data.x, aes_plot) +
                    geom_boxplot() +
                    ggtitle(str_glue("prefilter {comp_var}: n = {nrow(cur_row$data.x)}"))
                fil <- ggplot(cur_row$data.y, aes_plot) +
                    geom_boxplot() +
                    ggtitle(str_glue("postfilter {comp_var}: n = {nrow(cur_row$data.y)}"))
                return(unf + fil)
            })
            return(
                wrap_plots(var_plots, ncol = length(comp_vars)) +
                plot_annotation(title = str_glue("{unique(cur_row$data.x$ct_subcluster)} cluster_cell_type == cell_type"))
            )
        })

    pdf(file.path(out_path_base, "meta_comp.pdf"), height = 7, width = 7 * length(comp_vars))
    walk(comp_plots, print)
    dev.off()
}

if (!interactive()) {
    main()
}
