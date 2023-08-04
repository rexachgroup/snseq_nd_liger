set.seed(0)
liblist <- c("Seurat", "tidyverse", "variancePartition", "batchtools", "edgeR", "BiocParallel", "future.apply", "lme4", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
in_seurat_wk <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"
out_dir <- "../../analysis/seurat_lchen/seurat_cluster_lme/summary/"

main <- function() {
    dir.create(out_dir)
    cluster_wk <- readRDS(in_seurat_wk)
    cluster_wk <- cluster_wk %>% mutate(dge_filt = map(broom_join, filter_dx_specific_dge))
    cluster_wk %>%
        pluck("dge_filt") %>%
        bind_rows %>%
        write_csv(file.path(out_dir, "dge_filter.csv"))
    cluster_wk <- cluster_wk %>%
        mutate(dge_ct = map(broom_join, summarize_signif_dge))
    cluster_wk %>%
        pluck("dge_ct") %>%
        bind_rows %>%
        write_csv(file.path(out_dir, "dge_singif_count.csv"))

}

# split dx from broom_join results into separate column.
dge_extract_dx <- function(dge) {
    dge %>%
        pivot_longer(
            cols = -c("gene", "model", "region", "cell_type", "ct_cluster"), 
            names_pattern = "clinical_dx(.*)") %>%
        separate_wider_delim(
            cols = c("name"),
            delim = ".",
            names = c("dx", "dge_type"),
            too_many = "merge"
        ) %>%
        pivot_wider(
            names_from = c("dge_type"),
            values_from = c("value")
        )
}

summarize_signif_dge <- function(dge) { 
    dge_long <- dge_extract_dx(dge)
    dge_sum <- dge_long %>%
        group_by(ct_cluster, region, dx) %>%
        mutate(
            signif_up = estimate > 0.1 & p.value.adj < 0.1, 
            signif_dn = estimate < 0.1 & p.value.adj < 0.1
        ) %>%
        summarize(
            signif_up = sum(signif_up), 
            signif_dn = sum(signif_dn)
        )
    return(dge_sum)
}

filter_dx_specific_dge <- function(dge) {
    dge_long <- dge_extract_dx(dge)

    # genes that are significant in exactly one dx.
    dge_signif_filter <- dge_long %>%
        filter(p.value.adj < 0.1 & abs(estimate) > 0.1) %>%
        group_by(gene) %>%
        summarize(n = n()) %>%
        filter(n == 1)

    # genes that are not significant in at least two dx.
    dge_pval_nonsignif <- dge_long %>%
        filter(p.value > 0.05) %>%
        group_by(gene) %>%
        summarize(n = n()) %>%
        filter(n == 2)


    genes_dx_specific <- inner_join(dge_signif_filter, dge_pval_nonsignif, by = "gene")
    dge_f <- dge %>% filter(gene %in% genes_dx_specific$gene)
    return(dge_f)
}

main()
