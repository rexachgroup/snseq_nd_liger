# lchen
# aggregation of liger subclusters into one seurat object.

liblist <- c("tidyverse", "Seurat")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

LIGER_DIR <- "../analysis/liger_subcluster/20200824/"
OUT_DIR <- "../analysis/seurat_lchen/"

rdats <- list.files(LIGER_DIR, pattern = "*.rdat", full.names = TRUE)
rdats_filtered <- rdats[grep("test", rdats, invert = TRUE)]

merge_import <- function(rdats) {
    load(rdats[[1]])
    big_seurat <- liger_so
    rm(liger_so)

    for (r_path in rdats[2:length(rdats)]) {
        load(r_path)
        big_seurat <- merge(big_seurat, liger_so)
        rm(liger_so)
        gc()
    }
    return(big_seurat)
}

liger_so <- merge_import(rdats_filtered)

saveRDS(liger_so, file.path(OUT_DIR, "liger_subcluster_merged.rds"), compress = FALSE)
