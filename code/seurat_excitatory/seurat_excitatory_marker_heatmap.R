# Plot expression of liger clusters when stratified by excitatory layers.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "RColorBrewer", "scales")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"
options(future.globals.maxSize = Inf)

clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
excitatory_markers <- read_csv("../../resources/excitatory_markers_20191023.csv")
excitatory_markers2 <- read_csv("../../resources/excitatory_layers_20210318.csv")
out_dir <- "../../analysis/seurat_lchen/seurat_excitatory_layers/"
seurat_excitatory_meta <- readRDS("../../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds")
liger_meta <- readRDS("../../analysis/seurat_lchen/liger_subcluster_metadata.rds")
in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"

main <- function() {
    liger_assign <- select(liger_meta, liger_clusters, cell_ids)
    meta <- as_tibble(inner_join(seurat_excitatory_meta, liger_assign, by = "cell_ids"))
    marker_tb <- excitatory_markers %>%
        filter(!is.na(gene_symbol)) %>%
        filter(!duplicated(gene_symbol))
    sobj <- readRDS(in_seurat_rds)
    excludes <- read_xlsx(clusters_exclude_file)

    excitatory_meta <- meta %>%
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-")) %>%
        filter(cluster_cell_type == "excitatory") %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster)
    sobj <- subset(sobj, cells = excitatory_meta$cell_ids)
    gc()

    excitatory_meta %>%
        group_by(region) %>%
        group_walk(function(.x, .y) {
            out_path <- file.path(out_dir, str_glue("seurat_excitatory_heatmap_{unique(.x$region)}.pdf"))
            writeLines(out_path)
            heatmap_worker(sobj, .x, out_path)
        }, .keep = TRUE)
}

heatmap_worker <- function(sobj, excitatory_meta, out_path) {
    clust_summary_tb <- excitatory_meta %>%
        group_by(ct_subcluster) %>%
        group_nest %>%
        mutate(expr_mean = map(data, function(meta) {
            expr_m <- get_expr(sobj, meta, marker_tb)
            return(colMeans(expr_m))
        }))

    clust_expr_tb <- clust_summary_tb %>%
        select(ct_subcluster, expr_mean) %>%
        mutate(gene_name = map(expr_mean, names)) %>%
        unnest(c(gene_name, expr_mean))
    clust_expr_m <- pivot_matrix(clust_expr_tb, "ct_subcluster", "expr_mean", "gene_name") %>% scale(center = TRUE, scale = TRUE)
    clust_expr_scale <- colorRamp2(
        c(min(clust_expr_m), 0, max(clust_expr_m)),
        c(muted("blue"), "white", muted("red"))
    )

    col_annot_df <- excitatory_meta %>%
        count(ct_subcluster, clinical_dx) %>%
        group_by(ct_subcluster) %>%
        mutate(percent_cells = n / sum(n)) %>%
        ungroup()
    col_annot_barplot_counts <- excitatory_meta %>% count(ct_subcluster) %>% pull(n)
    col_annot_m <- pivot_matrix(col_annot_df, "clinical_dx", "percent_cells", "ct_subcluster")
    col_annot_obj <- HeatmapAnnotation(
        percent_dx = col_annot_m,
        number_of_cells = anno_barplot(col_annot_barplot_counts),
        col = list(
            percent_dx = colorRamp2(c(0, 1), c("white", muted("blue")))
        )
    )

    row_annotation_df <- marker_tb %>%
        filter(gene_symbol %in% rownames(clust_expr_m)) %>%
        column_to_rownames("gene_symbol") %>%
        select(marker_for)

    row_colors <- colorRampPalette(brewer.pal(9, "Blues"))(length(unique(row_annotation_df$marker_for)))
    names(row_colors) <- unique(row_annotation_df$marker_for)

    row_annotation_obj <- rowAnnotation(
        df = row_annotation_df,
        col = list(
            marker_for = row_colors
        )
    )

    heatmap_obj <- Heatmap(
        clust_expr_m,
        name = "normalized expression\nz-score",
        col = clust_expr_scale,
        cluster_rows = FALSE,
        top_annotation = col_annot_obj,
        right_annotation = row_annotation_obj,
        row_title = str_glue("Marker genes excitatory_layers_20191023.csv"),
        column_title = "Liger clusters",
        column_title_side = c("bottom"),
        border = TRUE
    )

    plot_width <- ncol(clust_expr_m) * 0.1 + 10
    plot_height <- plot_width

    pdf(out_path, width = plot_width, height = plot_height)
    print(heatmap_obj)
    graphics.off()
    gc()
    
}

get_expr <- function(sobj, meta, marker_tb) {
    genes <- marker_tb[marker_tb$gene_symbol %in% rownames(sobj), ]
    expr_m <- FetchData(sobj, vars = genes$gene_symbol, cells = meta$cell_ids)
    return(expr_m)
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(matches(paste0(c(cols_from, values_from, rows_from), collapse = "|"))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = values_from) %>%
        glimpse
    tb_matrix <- select(tb_pivot, -all_of(rows_from)) %>%
        as.matrix()
    rownames(tb_matrix) <- tb_pivot %>% rowwise() %>% summarize(paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>% pluck(1)
    glimpse(tb_matrix)
    return(tb_matrix)
}


