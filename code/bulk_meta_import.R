# Import and format bulk-seq data plus neurodegernation scores.
library(tidyverse)
in_bulk_meta <- "../resources/individual_nd.rds"
in_seurat_meta <- "../analysis/pci_import/pci_seurat_meta.rds"
out_dir <- "../analysis/bulk_meta/"

main <- function() { 
    dir.create(out_dir)
    bulk_meta <- readRDS(in_bulk_meta)
    seurat_meta <- readRDS(in_seurat_meta)
    bulk_meta_ref <- bulk_meta %>%
        mutate(
            Autopsy.ID = factor(paste0("P", Autopsy.ID)),
            type = factor(type),
            Region = fct_recode(Region, "calcarine" = "Calcarine", "insula" = "mid-insula", "preCG" = "PreCG")
        ) %>%
        group_by(Autopsy.ID, type, Region) %>%
        print

    seurat_ref <- seurat_meta %>%
        group_by(library_id, autopsy_id, region) %>%
        summarize(n = n(), .groups = "drop") %>%
        print() %>%
        select(-n)
    
    # match on autopsy id + region to get single cell library id labels
    nd_tb <- bulk_meta_ref %>%
        summarize(score = unique(score), .groups = "drop") %>%
        select(Autopsy.ID, type, score, Region) %>%
        rename(autopsy_id = Autopsy.ID, region = Region) %>%
        inner_join(seurat_ref, nd_tb, by = c("autopsy_id", "region"), multiple = "all") %>%
        print

    saveRDS(bulk_meta_ref, file.path(out_dir, "bulk_meta.rds"))
    saveRDS(nd_tb, file.path(out_dir, "nd_tb.rds"))
}

main()
