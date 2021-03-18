# Classify excitatory neurons based on layer markers. Plot + save resulting metadata.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "batchtools", "patchwork", "ggrastr", "ragg")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

excitatory_markers <- "../resources/excitatory_layers_20210318.csv"
in_seurat_rds <- "../analysis/pci_import/pci_seurat.rds"
out_dir <- "../analysis/seurat_lchen/seurat_excitatory_layers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pal_stallion = c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

main <- function() {
    sobj <- readRDS(in_seurat_rds)
    markers_tb <- read_csv(excitatory_markers) %>%
        filter(gene_symbol %in% rownames(sobj))

    sobj <- assign_layertype(sobj, markers_tb)
    sobj <- assign_cluster_layertype(sobj)
    sobj$cluster_celltype_layer <- ifelse(sobj$cluster_cell_type == "excitatory", sobj$cluster_layer_type, sobj$cluster_cell_type)
    sobj$celltype_layer <- ifelse(sobj$cell_type == "excitatory", sobj$layer_type, sobj$cell_type)
    sobj_umap <- as.data.frame(Embeddings(sobj, reduction = "umap"))

    sobj@meta.data %>% as_tibble %>%
        group_by(celltype_layer, cluster_ids) %>%
        summarize(n = n()) %>%
        print(n = Inf) %>%
        write_csv(file.path(out_dir, "ct_layer.csv"))
    pdf(file.path(out_dir, "excitatory_markers.pdf"), width = 10, height = 7)
    ggplot(sobj_umap, aes(x = UMAP_1, y = UMAP_2, color = celltype_layer)) + rasterize(geom_point(size = 1), dpi = 300, dev = "ragg_png") +
        scale_color_manual(values = setNames(pal_stallion, NULL)) +
        theme(aspect.ratio = 1)
    graphics.off()

    saveRDS(sobj@meta.data, file.path(out_dir, "sobj_celltype_meta.rds"), compress = FALSE)
}

assign_layertype <- function(sobj, markers, type_col, type_score) {
    writeLines("assign_celltype")
    markers <- group_by(markers, marker_group)
    ct_names <- group_keys(markers)$marker_group
    marker_split <- markers %>%
        group_split() %>%
        map(pluck("gene_symbol"))
    sobj@meta.data <- sobj@meta.data %>%
        dplyr::select(!contains(ct_names))
    sobj <- AddModuleScore(
        sobj,
        features = marker_split,
        assay = "RNA",
        name = paste0(ct_names, "fff")
    )

    # Rename columns by dropping trailing digit.
    sobj@meta.data <- sobj@meta.data %>%
        rename_with(.cols = contains(ct_names), ~str_extract(., "^.*(?=fff\\d+)"))

    max_ct_score <- sobj@meta.data %>%
        as_tibble %>%
        dplyr::select(cell_ids, contains(ct_names)) %>%
        pivot_longer(cols = contains(ct_names), names_to = "layer_type", values_to = "layer_score") %>%
        group_by(cell_ids) %>%
        dplyr::slice_max(layer_score, with_ties = FALSE) %>%
        dplyr::select(layer_type, layer_score)

    sobj@meta.data <- inner_join(sobj@meta.data, max_ct_score, by = "cell_ids")
    return(sobj)
}

assign_cluster_layertype <- function(sobj, cluster_ct_cutoff = 0.2) {
    writeLines("assign_cluster_layertype")

    cluster_counts <- sobj@meta.data %>%
        as_tibble %>%
        group_by(cluster_ids) %>%
        add_count(name = "layer_total") %>%
        group_by(cluster_ids, layer_total, layer_type) %>%
        count(name = "ct_total")

    cluster_assign <- cluster_counts %>%
        mutate(ct_frac = ct_total / layer_total) %>%
        group_by(cluster_ids) %>%
        slice_max(ct_frac, with_ties = FALSE) %>%
        mutate(cluster_layer_type = ifelse(ct_frac > cluster_ct_cutoff, layer_type, "unknown")) %>%
        dplyr::select(cluster_ids, cluster_layer_type)

    sobj@meta.data <- inner_join(sobj@meta.data, cluster_assign, by = "cluster_ids")
    return(sobj)
}

if (!interactive()) {
    main()
}
