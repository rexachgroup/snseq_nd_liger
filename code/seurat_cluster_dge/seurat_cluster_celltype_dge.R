# DGE of MAPT locus across celltype, per dx.

liblist <- c("Seurat", "tidyverse", "lme4", "lmerTest", "broom.mixed")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

plot_ts <- function() str_glue("seurat_cluster_celltype_dge.R {tools::md5sum('seurat_cluster_celltype_dge.R')} {date()}")

SEURAT <- "../../analysis/pci_import/pci_seurat.rds"
GENE_LIST <- "../../resources/cwow_genelist_WITHNOTES.csv"
OUT_DIR <- "../../analysis/seurat_lchen/seurat_cluster_celltype_dge/"

main <- function() {  
    dir.create(OUT_DIR)
    # Read genelist.
    marker_tb <- read_csv(GENE_LIST)

    sobj <- readRDS(SEURAT)
    subset_sobj <- subset(sobj, features = marker_tb$gene)
    rm(sobj)
    gc()
   
    marker_tb <- filter(marker_tb, gene %in% intersect(gene, rownames(subset_sobj)))

    dx_split <- subset_sobj@meta.data %>%
        as_tibble(rownames = "barcode") %>%
        group_by(clinical_dx) %>% group_nest()
    
    # subset to marker tb genes.
    dx_split <- dx_split %>%
        mutate(sobj_ct = map(data, ~subset(subset_sobj, cells = .$barcode)))

    # FindMarkers between excitatory and other cell types.
    dx_split <- dx_split %>%
        mutate(dge = map(sobj_ct, dx_ct_dge))

    # combine per-celltype results into one tb.
    dx_tb <- dx_split %>%
        mutate(dge_list = map(dge, function(dge_ct_list) {
            dge_ct_list %>%
                map(~rownames_to_column(.,"gene")) %>%
                bind_rows(.id = "ident.1") %>%
                mutate(ident.2 = "excitatory") %>%
                select(ident.1, ident.2, everything()) %>%
                rename(p_val_adj_celltype = p_val_adj)
        })) %>%
        select(-data, -sobj_ct, -dge) %>%
        unnest(dge_list)
    write_csv(dx_tb, file.path(OUT_DIR, "within_dx_celltype_dge.csv"))
}

dx_ct_dge <- function(sobj) {
    ident_base <- "excitatory"
    celltypes <- setdiff(unique(sobj$cluster_cell_type), ident_base)
    map(celltypes, function(ct) {
        FindMarkers(sobj,
            group.by = "cluster_cell_type",
            ident.1 = ct,
            ident.2 = ident_base,
            logfc.threshold = 0,
            min.pct = 0,
            test.use = "wilcox"
        )
    }) %>% setNames(nm = celltypes)
}

main()
