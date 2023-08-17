set.seed(0)
liblist <- c("tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

OUT_DIR <- "../../analysis/seurat_lchen/cellbender/dge_comp/"
SEURAT_DGE <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"
CELLBENDER_DGE <- "../../analysis/seurat_lchen/cellbender/dge/cluster_lme.rds"
SEURAT_SUMMARY <- "../../analysis/seurat_lchen/seurat_cluster_lme/dx_dge_summary.csv"
CELLBENDER_SUMMARY <- "../../analysis/seurat_lchen/cellbender/dge/dx_dge_summary.csv"

# 1. read Seurat dge
# 2. get dx-specific gene esimate + region
# 3. join agianst cellbender dge results
# 4. diff estimates
main <- function() {
    dir.create(OUT_DIR)
    seurat_dge_tb <- readRDS(SEURAT_DGE)   
    cellbender_dge_tb <- readRDS(CELLBENDER_DGE)

    seurat_dge_tb <- seurat_dge_tb %>% 
        mutate(dge_filt = map(broom_join, filter_dx_specific_dge))

    seurat_filt_tb <- seurat_dge_tb %>%
        pluck("dge_filt") %>%
        bind_rows

    seurat_filt_long_tb <- dge_extract_dx(seurat_filt_tb)

    cellbender_long_tb <- cellbender_dge_tb %>%
        dge_extract_dx


    s_join_tb <- seurat_filt_long_tb %>%
        group_by(gene, ct_cluster) %>%
        group_nest()

    c_join_tb <- cellbender_long_tb %>%
        group_by(gene, ct_cluster) %>%
        group_nest()

    dge_join_tb <- inner_join(s_join_tb, c_join_tb, by = c("gene", "ct_cluster"))

    dge_join_tb <- dge_join_tb %>%
        mutate(dge_diff = pmap(., function(...) {
            cr <- list(...)
            data_join <- left_join(cr$data.x, cr$data.y, by= c("dx"))
            data_join <- data_join %>%
                mutate(estimate_diff = estimate.y - estimate.x) %>%
                select(model = model.x, region = region.x, dx, estimate_diff)
        }))
    dge_join_tb %>%
        select(gene, ct_cluster, dge_diff) %>%
        unnest(dge_diff) %>%
        write_csv(file.path(OUT_DIR, "dge_estimate_diff.csv")) 
}

dge_extract_dx <- function(dge) {
    dge %>%
        pivot_longer(
            cols = -c("gene", "model", "region", "cell_type", "ct_cluster"), 
            names_pattern = "clinical_dx(.*)") %>%
        separate(
            col = c("name"),
            sep = "\\.",
            into = c("dx", "dge_type"),
            extra = "merge"
        ) %>%
        pivot_wider(
            names_from = c("dge_type"),
            values_from = c("value")
        )
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

if (!interactive()) main()
