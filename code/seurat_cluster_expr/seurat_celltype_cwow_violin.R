# Calculate pseudobulk per sample for cwow genelist.
# Plot as violin plot.

liblist <- c("Seurat", "ggplot2", "tidyverse", "patchwork", "ragg")
l <- lapply(liblist, library, character.only = TRUE, quietly = TRUE)

GENE_LIST <- "../../resources/cwow_genelist_WITHNOTES.csv"
SEURAT <- "../../analysis/pci_import/pci_seurat.rds"
OUT_DIR <- "../../analysis/seurat_lchen/seurat_cluster_expr/"

plot_ts <- function() str_glue("seurat_celltype_cwow_violin.R {tools::md5sum('seurat_celltype_cwow_violin.R')} {date()}")

main <- function() {
    dir.create(OUT_DIR, recursive = T)
    marker_tb <- read_csv(GENE_LIST)

    sobj <- readRDS(SEURAT)
    subset_sobj <- subset(sobj, features = marker_tb$gene)
    rm(sobj)
    gc()

    marker_tb <- filter(marker_tb, gene %in% rownames(subset_sobj))

    celltype_split <- subset_sobj@meta.data %>%
        as_tibble(rownames = "barcode") %>%
        group_by(cluster_cell_type) %>% group_nest()

    celltype_split <- celltype_split %>%
        mutate(sobj_ct = map(data, ~subset(subset_sobj, cells = .$barcode)))

    celltype_split <- celltype_split %>%
        mutate(sobj_ct_test = map(sobj_ct, dx_dge))
    
    celltype_split <- celltype_split %>%
        mutate(pseudo_bulk = map(sobj_ct, pseudo_bulk_subj))

    celltype_split <- celltype_split %>%
        mutate(violin_plotlist = pmap(., function(...) {
            cr <- list(...)
            violin_genelist_patchwork(cr$sobj_ct, cr$pseudo_bulk, marker_tb$gene)
        }))

    pwalk(celltype_split, function(...) {
        cr <- list(...)
        ndim <- ceiling(sqrt(length(marker_tb$gene))) * 3
        out_path <- str_glue("{OUT_DIR}/{cr$cluster_cell_type}.pdf")
        writeLines(out_path)
        pdf(out_path, width = ndim, height = ndim)
        print(cr$violin_plotlist + 
              plot_annotation(title = str_glue("{cr$cluster_cell_type} average expression per library"), subtitle = plot_ts())
        )
        graphics.off()
    })
    celltype_split %>%
        select(-data, -sobj_ct, -pseudo_bulk) %>%
        unnest(sobj_ct_test) %>%
        write_csv(file.path(OUT_DIR, "seurat_dx_test.csv"))
}

# dge per dx v. control.
# p_val_adj is based on total gene list of cluster before filtering (n = 33)
dx_dge <- function(sobj) {
    ident_col = "clinical_dx"
    ident_base = "Control"
    test_idents <- setdiff(unique(sobj[[ident_col, drop = T]]), ident_base)
    map(test_idents, function(ident_test) {
        FindMarkers(sobj,
            group.by = ident_col,
            ident.1 = ident_test,
            ident.2 = ident_base,
            logfc.threshold = 0,
            min.pct = 0,
            test.use = "wilcox",
        ) %>% as_tibble(rownames = "gene")
    }) %>% setNames(nm = test_idents) %>% bind_rows(.id = "ident.1") %>%
    mutate(ident.2 = ident_base)
}
violin_genelist_patchwork <- function(sobj, pseudo_bulk, genes) {
    genes <- genes[genes %in% rownames(sobj)]
    ndim <- ceiling(sqrt(length(genes)))
    plotlist <- map(genes, ~violin_gene_subj_dx(sobj, pseudo_bulk, .))
    return(wrap_plots(plotlist, ncol = ndim, nrow = ndim))
}

violin_gene_dx <- function(sobj, gene) {
    return(VlnPlot(sobj, features = gene, group.by = "clinical_dx", fill.by = "clinical_dx"))
}

violin_gene_subj_dx <- function(sobj, pseudo_bulk, gene) {
    plot_tb <- sobj@meta.data %>%
        as_tibble(rownames = "barcode") %>%
        group_by(library_id, clinical_dx) %>%
        summarize(n = n(), .groups = "drop") %>%
        mutate(aveExpr = pseudo_bulk$RNA[gene, library_id])

    gg <- ggplot(plot_tb, aes(x = clinical_dx, y = aveExpr, fill = clinical_dx)) +
        geom_violin() +
        geom_jitter(width = 0.3, height = 0) +
        labs(title = gene) +
        theme(legend.position = "none")

    return(gg)
}

pseudo_bulk_subj <- function(sobj) {
    return(AverageExpression(sobj, group.by = "library_id"))
}

if (!interactive()) main()
