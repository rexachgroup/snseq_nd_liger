set.seed(0)
liblist <- c("tidyverse", "patchwork", "ggrastr")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

OUT_DIR <- "../../analysis/seurat_lchen/cellbender/dge_comp_no_logumi/"
SEURAT_DGE <- "../../analysis/seurat_lchen/seurat_cluster_lme/cluster_wk.rds"
CELLBENDER_DGE <- "../../analysis/seurat_lchen/cellbender/dge_no_logumi/cluster_lme.rds"
SEURAT_SUMMARY <- "../../analysis/seurat_lchen/seurat_cluster_lme/dx_dge_summary.csv"
CELLBENDER_SUMMARY <- "../../analysis/seurat_lchen/cellbender/dge_no_logumi/dx_dge_summary.csv"

main <- function() {
    dir.create(OUT_DIR)
    seurat_dge_tb <- readRDS(SEURAT_DGE)   
    cellbender_dge_tb <- readRDS(CELLBENDER_DGE)

    seurat_dge_long_tb <- seurat_dge_tb %>%
        mutate(dge_filt = map(broom_join, dge_extract_dx)) %>%
        pluck("dge_filt") %>%
        bind_rows

    cellbender_dge_long_tb <- dge_extract_dx(cellbender_dge_tb)
 
    dx_levels <- unique(seurat_dge_long_tb$dx)

    dx_specific_seurat <- map(dx_levels, ~mark_dx_specific_dge(seurat_dge_long_tb, .))
    dx_specific_cellbender <- map(dx_levels, ~mark_dx_specific_dge(cellbender_dge_long_tb, .))

    dx_specific_seurat_long <- bind_rows(dx_specific_seurat)
    dx_specific_cellbender_long <- bind_rows(dx_specific_cellbender)

    # join seurat / cellbender dge
    dge_join_tb <- full_join(seurat_dge_long_tb, cellbender_dge_long_tb, by = c("gene", "ct_cluster", "dx"), suffix = c(".seurat", ".cellbender"))

    dge_diff_tb <- dge_join_tb %>%
        filter(!is.null(estimate.seurat) & !is.null(estimate.cellbender))
#   dge_join_tb <- dge_join_tb %>%
#       mutate(dge_diff = pmap(., function(...) {
#           cr <- list(...)
#           data_join <- left_join(
#               cr$data.x, cr$data.y, 
#               by = c("dx"),
#               suffix = c(".seurat", ".cellbender")
#           )
#           data_join <- data_join %>%
#               mutate(estimate_diff = estimate.cellbender - estimate.seurat)
#       }))

#   # write out dge table
#   dge_diff_tb <- dge_join_tb %>%
#       select(gene, ct_cluster, dge_diff) %>%
#       unnest(dge_diff)
#   write_csv(dge_diff_tb, file.path(OUT_DIR, "dge_estimate_diff.csv")) 

    dx_specific_seurat_nest <- dx_specific_seurat_long %>%
        group_by(gene, ct_cluster, dx) %>%
        group_nest(.key = "seurat_specific")
    dx_specific_cellbender_nest <- dx_specific_cellbender_long %>%
        group_by(gene, ct_cluster, dx) %>%
        group_nest(.key = "cellbender_specific")

    dge_dx_comp_tb <- dge_join_tb %>%
        left_join(dx_specific_seurat_nest) %>%
        left_join(dx_specific_cellbender_nest)

    saveRDS(dge_dx_comp_tb, file.path(OUT_DIR, "dge_estimate_cache.rds"), compress = F)

    dge_dx_comp_tb <- dge_dx_comp_tb %>%
        mutate(dge_color_data = unlist(pmap(., function(...) {
            cr <- list(...)
            seurat_specific <- !is.null(cr$seurat_specific)
            cellbender_specific <- !is.null(cr$cellbender_specific)
            if (seurat_specific & cellbender_specific) {
                signif_code <- "both_specific"
            } else if (seurat_specific) {
                signif_code <- "seurat_specific"
            } else if (cellbender_specific) {
                signif_code <- "cellbender_specific"
            } else {
                signif_code <- "none"
            }
            return(signif_code)
        })))
    
    dge_dx_comp_tb <- dge_dx_comp_tb %>%
        mutate(dge_alpha_data = unlist(pmap(., function(...) {
            cr <- list(...)
            seurat_specific <- !is.null(cr$seurat_specific)
            cellbender_specific <- !is.null(cr$cellbender_specific)
            if (seurat_specific && cellbender_specific) {
                signif_code <- 1
            } else if (seurat_specific) {
                signif_code <- 1
            } else if (cellbender_specific) {
                signif_code <- 1
            } else {
                signif_code <- 0.1
            }
            return(signif_code)
        })))

    dge_dx_plot_tb <- dge_dx_comp_tb %>%
        mutate(dge_color_data = fct_relevel(dge_color_data, c("none", "seurat_specific", "cellbender_specific", "both_specific"))) %>%
        replace_na(list(
            statistic.seurat = 0,
            statistic.cellbender = 0
        )) %>%
        group_by(ct_cluster, dx) %>%
        group_nest()

    colormap <- setNames(scales::brewer_pal(palette = "Set1")(4), 
        c("none", "seurat_specific", "cellbender_specific", "both_specific"))

    write_csv(dge_dx_comp_tb, file.path(OUT_DIR, "dge_dx_comp_tb.csv"))
    pdf(file.path(OUT_DIR, "dge_statistic_scatter.pdf"))
    dge_dx_plot_tb %>%
        pwalk(., function(...) {
            cr <- list(...)
            print(cr$data)
            dat <- cr$data %>%
                arrange(dge_color_data)
            gg <- ggplot(dat) +
                ggrastr::rasterize(list(
                    geom_point(mapping = aes(x = statistic.seurat, y = statistic.cellbender, color = dge_color_data, alpha = dge_alpha_data), size = 0.5),
                    scale_color_manual(values = colormap, breaks = names(colormap), drop = F),
                    scale_alpha()
                ), dpi = 100) + 
                labs(title = str_glue("statistic comparison, {cr$ct_cluster} {cr$dx}"))
            print(gg)
        })
    graphics.off()
    
    pdf(file.path(OUT_DIR, "dge_estimate_scatter.pdf"), width = 10, height = 10)
    dge_dx_plot_tb %>%
        pwalk(., function(...) {
            cr <- list(...)
            print(cr$data)
            dat <- cr$data %>%
                arrange(desc(dge_color_data))
            gg <- ggplot(dat, aes(x = estimate.seurat, y = estimate.cellbender, color = dge_color_data, alpha = dge_alpha_data)) +
                ggrastr::rasterize(list(
                    geom_point(size = 0.5),
                    scale_color_manual(values = colormap)
                ), dpi = 100) + 
                labs(title = str_glue("estimate comparison, {cr$ct_cluster} {cr$dx}"))
            print(gg)
            
        })
    graphics.off() 
}


comp_plot <- function(dge_tb) {
    ggplot(dge_tb) +
        geom_point(aes(x = statistic.seurat, y = statistic.cellbender)) +
        facet_wrap("ct_cluster")


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

mark_dx_specific_dge <- function(dge_long, tgt_dx) {
    # genes that are significant in exactly one dx.
    dge_signif_filter <- dge_long %>%
        filter(p.value.adj < 0.1 & abs(estimate) > 0.1 & abs(statistic) > 2 & dx == tgt_dx) %>%
        group_by(gene, ct_cluster) %>%
        summarize(n = n(), dx = unique(dx), .groups = "keep") %>%
        filter(n == 1)

    # genes that are not significant in at least two dx.
    dge_pval_nonsignif <- dge_long %>%
        filter(p.value > 0.05 & dx != tgt_dx) %>%
        group_by(gene, ct_cluster) %>%
        summarize(n = n(), .groups = "keep") %>%
        filter(n == 2)

    genes_dx_specific <- inner_join(dge_signif_filter, dge_pval_nonsignif, by = c("gene", "ct_cluster"))
    return(genes_dx_specific)
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

filter_dx_specific_dge_long <- function(dge) {

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
