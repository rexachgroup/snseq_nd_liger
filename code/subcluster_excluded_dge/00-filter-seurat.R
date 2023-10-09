set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

in_liger_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_seurat_sobj <- "../../analysis/pci_import/pci_seurat.rds"
in_cellbender_sobj <- "../../analysis/"
cellbender_tb <- "../../analysis/seurat_lchen/cellbender/meta_bc_join.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_excluded_dge/seurat_filter/"

main <- function() {
    dir.create(OUT_DIR, recursive = T)
    sobj <- readRDS(in_seurat_sobj)   
    liger_meta <- readRDS(in_liger_meta)
    excludes <- read_xlsx(clusters_exclude_file)
    cellbender_bc <- readRDS(cellbender_tb)

    cluster_celltype_v_celltype_breakdown(cellbender_bc, file.path(OUT_DIR, "00_unfiltered"))

    excludes_f <- cellbender_bc %>% filter(
            !ct_subcluster %in% excludes$ct_subcluster, 
            cluster_cell_type == cell_type,
            !cell_type %in% c("t_cell", "ependymal")
        )

    cluster_celltype_v_celltype_breakdown(excludes_f, file.path(OUT_DIR, "01_subcluster_filtered")) 

    excludes_cellbender <- cellbender_bc %>%
        filter(
            !ct_subcluster %in% excludes$ct_subcluster,
            cluster_cell_type == cell_type
            latent_cell_probability >= 0.5,
            !cell_type %in% c("t_cell", "ependymal")
        )
        
    cluster_celltype_v_celltype_breakdown(excludes_cellbender, file.path(OUT_DIR, "02_subcluster_cellbender_filtered"))

    tmp <- cellbender_bc %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        filter(cluster_cell_type == "excitatory") %>%
        mutate(matching_cluster_cell_type = cluster_cell_type == cell_type)
    pdf(file.path(OUT_DIR, "cluster_celltype_dist.pdf"))
    ggplot(tmp, aes(x = cell_type, y = number_umi)) +
        geom_boxplot()
    graphics.off()

}

cluster_celltype_v_celltype_breakdown <- function(sobj_meta, out_prefix, filter_statement) {
   sobj_meta %>%
       select(cluster_cell_type, cell_type) %>%
       table() %>%
       as.matrix() %>%
       write.csv(str_glue("{out_prefix}_ct_matrix.csv"))
   sobj_meta %>%
       group_by(cluster_cell_type) %>%
        summarize(frac_matching = sum(cell_type == cluster_cell_type) / n()) %>%
        write.csv(str_glue("{out_prefix}_ct_table.csv"))
    cell_types <- unique(sobj_meta$cluster_cell_type)
    map(cell_types, function(x) {
        meta_c <- sobj_meta %>%
            filter(cluster_cell_type == x) %>%
            mutate(
                matching_cluster_cell_type = cluster_cell_type == cell_type,
                removed_by_cellbender = latent_cell_probability < 0.5
            )
        meta_c %>%
            select(matching_cluster_cell_type, removed_by_cellbender) %>% 
            table %>% 
            as.data.frame %>%
            mutate(cluster_cell_type = x) %>%
            select(cluster_cell_type, everything())
    }) %>%
    bind_rows %>%
    write.csv(str_glue("{out_prefix}_ct_v_cellbender.csv"))
    write.csv(sobj_meta, str_glue("{out_prefix}_meta.csv"))
}

