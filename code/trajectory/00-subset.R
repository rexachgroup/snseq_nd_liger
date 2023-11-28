# Subset to insula layer 2/3 IT for trajectory analysis.
set.seed(0)
options(deparse.max.lines = 5)
liblist <- c("Seurat", "tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

SEURAT_OBJ <- "../../analysis/pci_import/pci_seurat.rds"
META <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
#AZIMUTH_META <- "../../ext/azimuth/data/meta.csv"
OUT_DIR <- "../../analysis/seurat_lchen/trajectory/insula-subset/"

SUBCLUSTER_LIST <- c(
    "insula-excitatory-2",
    "insula-excitatory-5",
    "insula-excitatory-6"
)

main <- function() {
    dir.create(OUT_DIR)
    sobj <- readRDS(SEURAT_OBJ)
    meta <- readRDS(META)

    meta_subset <- meta %>%
        filter(region == "insula", ct_subcluster %in% SUBCLUSTER_LIST, cell_ids %in% colnames(sobj))
    meta_df <- meta_subset %>%
        as.data.frame
    rownames(meta_df) <- meta_df$cell_ids

    subset_sobj <- subset(sobj, cells = meta_subset$cell_ids)
    subset_sobj <- AddMetaData(subset_sobj, meta_df)
    saveRDS(subset_sobj, file.path(OUT_DIR, "sobj.rds"), compress = F)
}

if (!interactive()) main()
