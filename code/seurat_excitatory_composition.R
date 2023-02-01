# Per-dx celltype proportions.
liblist <- c("tidyverse", "Seurat", "ggpubr", "broom", "rstatix")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

in_celltype_meta <- readRDS("../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds")
out_dir <- "../analysis/seurat_lchen/seurat_excitatory_layers"

plot_ts <- function() str_glue("{date()} seurat_excitatory_composition.R {tools::md5sum('seurat_excitatory_composition.R')}")

main <- function() {

    # Group based on cluster_celltype_layer.
    cluster_layer_dat <- in_celltype_meta %>%
        group_by(clinical_dx, region) %>%
        ungroup %>%
        mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
        group_by(library_id, autopsy_id, clinical_dx, region, cluster_celltype_layer) %>%
        summarise(n = n(), percent_mito = mean(percent_mito), age = first(age), pmi = first(pmi), sex = first(sex), .groups = "drop") %>%
        group_by(library_id) %>%
        mutate(percent_of_library = n / sum(n)) %>%
        arrange(region, clinical_dx, cluster_celltype_layer)

    # dx and region-specific summaries.
    region_sum <- cluster_layer_dat %>%
        group_by(region) %>%
        summarize(n = sum(n))

    dx_sum <- cluster_layer_dat %>%
        group_by(clinical_dx) %>%
        summarize(n = sum(n))
    
    dxreg_sum <- cluster_layer_dat %>%
        group_by(clinical_dx, region) %>%
        summarize(n = sum(n))
    
    aut_sum <- cluster_layer_dat %>%
        group_by(autopsy_id) %>%
        summarize(n = sum(n))
    
    lib_sum <- cluster_layer_dat %>%
        group_by(library_id) %>%
        summarize(n = sum(n))
    
    celltype_sum <- cluster_layer_dat %>%
        group_by(cluster_celltype_layer) %>%
        summarize(n = sum(n))

    write.xlsx(
        list(
            full_table = cluster_layer_dat,
            group_by_celltype = celltype_sum,
            group_by_region = region_sum,
            group_by_dx = dx_sum,
            group_by_dx_region = dxreg_sum,
            group_by_autopsy_id = aut_sum,
            group_by_library_id = lib_sum,
            ts = plot_ts()
        ),
        file.path(out_dir, "cluster_celltype_layer_proportions.xlsx")
    )
}

if (!interactive()) {
    main()
}
