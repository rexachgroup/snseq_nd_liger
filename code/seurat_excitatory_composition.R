# Per-dx celltype proportions.
liblist <- c("tidyverse", "Seurat", "openxlsx")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)

in_celltype_meta <- "../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds"
in_nd_tb <- "../analysis/bulk_meta/nd_tb.rds"
out_dir <- "../analysis/seurat_lchen/seurat_excitatory_layers"

plot_ts <- function() str_glue("{date()} seurat_excitatory_composition.R {tools::md5sum('seurat_excitatory_composition.R')}")

main <- function() {
    celltype_meta <- readRDS(in_celltype_meta)
    nd_tb <- readRDS(in_nd_tb)
    nd_tb_wd <- nd_tb %>%
        pivot_wider(names_from = "type", values_from = "score")

    # Group based on cluster_celltype_layer.
    cluster_layer_dat <- celltype_meta %>%
        mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
        group_by(library_id, autopsy_id, clinical_dx, region, cluster_celltype_layer) %>%
        summarise(n = n(), percent_mito = mean(percent_mito), age = first(age), pmi = first(pmi), sex = first(sex), .groups = "drop") %>%
        group_by(library_id) %>%
        mutate(percent_of_library = n / sum(n)) %>%
        arrange(region, clinical_dx, cluster_celltype_layer) %>%
        left_join(nd_tb_wd, by = c("autopsy_id", "region"))

    # library summary.
    aut_sum <- cluster_layer_dat %>%
        group_by(autopsy_id, region) %>%
        summarize(n = sum(n),clinical_dx = first(clinical_dx), age = first(age), pmi = first(pmi), sex = first(sex), composite = first(composite), nd = first(nd), tau = first(tau))

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
    
    # library and celltype summaries
    lib_sum <- cluster_layer_dat %>%
        group_by(library_id) %>%
        summarize(n = sum(n))
    
    celltype_sum <- cluster_layer_dat %>%
        group_by(cluster_celltype_layer) %>%
        summarize(n = sum(n))

    write.xlsx(
        list(
            full_table = cluster_layer_dat,
            group_by_autopsy_id = aut_sum,
            group_by_celltype = celltype_sum,
            group_by_region = region_sum,
            group_by_dx = dx_sum,
            group_by_dx_region = dxreg_sum,
            group_by_library_id = lib_sum,
            ts = plot_ts()
        ),
        file.path(out_dir, "cluster_celltype_layer_proportions.xlsx")
    )
}

if (!interactive()) {
    main()
}
