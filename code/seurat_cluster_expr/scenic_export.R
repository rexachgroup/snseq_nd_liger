# export metadata for specific clusters identified in scenic
set.seed(0)
liblist <- c("Seurat", "tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
seurat_meta <- "../../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds"
bulkseq_meta <- "../../resources/individual_nd.rds"
scenic_insula <- "../../analysis/scenic/insula_nucseq-nd_seurat_object.rds"
scenic_precg <- "../../analysis/scenic/precg_nucseq-nd_seurat_object.rds"
out_dir <- "../../analysis/seurat_lchen/scenic_export/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

main <- function() {
    meta <- readRDS(seurat_meta)
    seurat_insula <- readRDS(scenic_insula)
    seurat_insula_meta <- seurat_insula@meta.data

    seurat_precg <- readRDS(scenic_precg)
    seurat_precg_meta <- seurat_precg@meta.data

    cid <- seurat_insula_meta %>%
        filter(cluster_ids == 31) %>%
        pluck("cell_ids")

    insula_counts <- inner_join(
        seurat_insula_meta %>%
            filter(cluster_ids %in% c(36)) %>%
            group_by(cluster_ids, library_id) %>%
            summarize(n = n()),
        seurat_insula_meta %>%
            filter(cluster_ids %in% c(36)) %>%
            group_by(cluster_ids, library_id) %>%
            slice_head(),
        by = c("cluster_ids", "library_id")
    )

    precg_counts <- inner_join(
        seurat_precg_meta %>%
            filter(cluster_ids %in% c(31, 24)) %>%
            group_by(cluster_ids, library_id) %>%
            summarize(n = n()),
        seurat_precg_meta %>%
            filter(cluster_ids %in% c(31, 24)) %>%
            group_by(cluster_ids, library_id) %>%
            slice_head(),
        by = c("cluster_ids", "library_id") 
    )

    write_csv(insula_counts, file.path(out_dir, "insula_scenic_counts.csv"))
    write_csv(precg_counts, file.path(out_dir, "precg_scenic_counts.csv"))
}
main()
