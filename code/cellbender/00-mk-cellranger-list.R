library(tidyverse)

AGGR_FILE <- "/u/project/geschwind/drewse/g_singlecell/cellranger/data/20200710/aggregate/outs/aggregation.csv"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/"
dir.create(OUT_DIR)

#filelist <- list.files(LS_DIR, recursive = T, full.names = T, pattern = "raw_feature_bc_matrix.h5")
#is_cellranger_aggr <- file.exists(file.path(dirname(filelist), "aggregation.csv"))

file_tb <- read_csv(AGGR_FILE)

file_tb <- file_tb %>%
    mutate(path_components = map(molecule_h5, ~unlist(str_split(., "/"))))
file_tb <- file_tb %>%
    mutate(
        len = unlist(map(path_components, ~length(unlist(.x))))
    )
file_tb <- file_tb %>%
    mutate(
        date = unlist(pmap(list(path_components, len), ~.x[[.y - 3]])),
        sample_id = unlist(pmap(list(path_components, len), ~.x[[.y - 2]])),
    )
write_csv(file_tb, file.path(OUT_DIR, "cellranger_h5.csv"))
