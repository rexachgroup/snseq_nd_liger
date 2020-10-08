# lchen
# dump metadata from each seurat object.

liblist <- c("tidyverse", "Seurat")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

LIGER_DIR <- "../analysis/liger_subcluster/20200824/"
OUT_DIR <- "../analysis/seurat_lchen/"

rdats <- list.files(LIGER_DIR, pattern = "*.rdat", full.names = TRUE)
rdats_filtered <- rdats[!grepl("test", rdats)]

liger_metadata <- map(rdats_filtered, function(rdat_path) {
    load(rdat_path)
    meta <- liger_so@meta.data %>%
        as_tibble(rownames = "UMI")
    meta$filepath <- rdat_path
    return(meta)
})

liger_tb <- bind_rows(liger_metadata)

saveRDS(liger_tb, file.path(OUT_DIR, "liger_subcluster_metadata.rds"), compress = FALSE)
write_csv(liger_tb, file.path(OUT_DIR, "liger_subcluster_meta.csv"))
gc()
