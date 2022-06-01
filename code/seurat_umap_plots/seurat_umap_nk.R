# Plot metadata columns overlaid on UMAP embedding.
# Plot individual nk gene names overlaid on UMAP embedding.

set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "patchwork", "ggrastr")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

out_dir <- "../../analysis/seurat_lchen/seurat_umap"
in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
in_seurat_meta <- "../../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds"
nk_markers <- "../../resources/Immune_markers_toplot_2021.xlsx"
sobj_umap_path <- file.path(out_dir, "sobj_meta_umap.rds")
sobj_umap_expr_subset_path <- file.path(out_dir, "sobj_expr_subset.rds")

plot_ts <- str_glue("{date()} {system('md5sum seurat_umap.R', intern = TRUE)}")
pal_stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

main <- function() {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(sobj_umap_path)) { extract_umap(in_seurat_rds, in_seurat_meta, sobj_umap_path) }
    if (!file.exists(sobj_umap_expr_subset_path)) { extract_nk_genes(in_seurat_rds, nk_markers, sobj_umap_expr_subset_path) }

    umap_meta <- readRDS(sobj_umap_path)
    nk_expr_mat <- readRDS(sobj_umap_expr_subset_path)

    pdf(file.path(out_dir, "sobj_umap.pdf"), width = 10, height = 7)
    print(ggplot_umap_dx(umap_meta))
    print(ggplot_umap_region(umap_meta))
    print(ggplot_umap_cluster(umap_meta))
    graphics.off()

    nk_list <- ggplot_umap_expr_mat(umap_meta, nk_expr_mat)
    plot_nrow <- ceiling(sqrt(length(nk_list)))
    plot_ncol <- ceiling(sqrt(length(nk_list)))
    pdf(file.path(out_dir, "sobj_nk_plots.pdf"), width = 5 * plot_nrow, height = 5 * plot_ncol)
    print(wrap_plots(nk_list, nrow = plot_nrow, ncol = plot_ncol) + plot_annotation(title = str_glue("{plot_ts} gene umap")))
    graphics.off()
}

extract_umap <- function(in_seurat_rds, in_seurat_meta, sobj_umap_path) {
    sobj <- readRDS(in_seurat_rds)
    meta <- readRDS(in_seurat_meta)
    sobj_umap <- Embeddings(sobj, "umap") %>%
        as.data.frame %>%
        rownames_to_column("cell_ids") %>%
        inner_join(meta, ., by = "cell_ids") 

    saveRDS(sobj_umap, sobj_umap_path, compress = FALSE)    
}

extract_nk_genes <- function(in_seurat_rds, nk_markers, sobj_umap_expr_subset_path) {
    sobj <- readRDS(in_seurat_rds)
    nk_marker_tb <- read_xlsx(nk_markers)
    nk_marker_list <- pluck(nk_marker_tb, "Genesymbol") %>% subset(., . %in% rownames(sobj))
    expr_mat <- GetAssayData(sobj, "data")[nk_marker_list, ]
    saveRDS(expr_mat, sobj_umap_expr_subset_path, compress = FALSE)
}

ggplot_umap_dx <- function(meta) {
    ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = clinical_dx)) +
        rasterize(geom_point(size = 1), dpi = 300, dev = "ragg_png") +
        scale_color_manual(values = setNames(pal_stallion, NULL)) +
        theme(aspect.ratio = 1) +
        ggtitle(str_glue("{plot_ts} clinical_dx umap"))
}

ggplot_umap_region <- function(meta) {
    ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = region)) +
        rasterize(geom_point(size = 1), dpi = 300, dev = "ragg_png") +
        scale_color_manual(values = setNames(pal_stallion, NULL)) +
        theme(aspect.ratio = 1) +
        ggtitle(str_glue("{plot_ts} region umap"))
}

ggplot_umap_cluster <- function(meta) {
    cluster_centers <- meta %>%
        group_by(cluster_celltype_layer) %>%
        summarize(cluster_x = mean(UMAP_1), cluster_y = mean(UMAP_2))
    ggplot(meta, aes(x = UMAP_1, y = UMAP_2, color = cluster_celltype_layer)) + 
        rasterize(geom_point(size = 0), dpi = 300, dev = "ragg_png") +
        geom_text(data = cluster_centers, mapping = aes(x = cluster_x, y = cluster_y, label = cluster_celltype_layer), colour = "black", fill = "black") +
        scale_color_manual(values = setNames(pal_stallion, NULL)) +
        theme(aspect.ratio = 1) +
        ggtitle(str_glue("{plot_ts} cluster umap"))
}

ggplot_umap_expr_mat <- function(meta, expr_mat) {
    map(rownames(expr_mat), function(gene_name) {
        expr_tb <- expr_mat[gene_name, ] %>%
            as.data.frame %>%
            rownames_to_column("cell_ids") %>%
            as_tibble %>%
            rename("gene_expr" = ".")

        plot_tb <- inner_join(meta, expr_tb, by = "cell_ids") %>%
            select(cell_ids, UMAP_1, UMAP_2, gene_expr) %>%
            mutate(alpha = ifelse(gene_expr > 0, 0.8, 0.2))
        gg <- ggplot(plot_tb, aes(x = UMAP_1, y = UMAP_2, color = gene_expr, alpha = alpha)) +
            rasterize(geom_point(size = 0), dpi = 200, dev = "ragg_png") +
            scale_color_distiller(type = "div", palette = "RdGy", direction = -1, na.value = "grey90") + 
            theme(aspect.ratio = 1) +
            ggtitle(str_glue("{gene_name}"))

        return(gg)
    })
}

if (!interactive()) {
    main()
}
