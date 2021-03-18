# Per-dx celltype proportions for all control samples.
# Downsampling per dx / region to get consistent propoertion estimates(?) between different regions.
liblist <- c("tidyverse", "Seurat", "ggpubr", "broom")
lapply(liblist, require, character.only = TRUE, quiet = TRUE)
commit_id <- str_glue("{system('git show --format=\"%h\" -s', intern = TRUE)}")

#in_seurat_metadata <- readRDS("../analysis/seurat_lchen/liger_subcluster_subset.rds")
in_celltype_meta <- readRDS("../analysis/seurat_lchen/seurat_excitatory_layers/sobj_celltype_meta.rds")
out_dir <- "../analysis/seurat_lchen/seurat_excitatory_layers"

main <- function() {
    set.seed(0)
    # Individual layer assignments.
    layer_dat <- in_celltype_meta %>%
        group_by(clinical_dx, region) %>%
        slice_sample(n = 20000) %>%
        filter(!celltype_layer %in% c("pericyte", "ependymal", "t_cell")) %>%
        ungroup %>%
        mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
        group_by(library_id, clinical_dx, region, celltype_layer) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(library_id) %>%
        mutate(percent_of_library = n / sum(n)) %>%
        arrange(region, clinical_dx, celltype_layer)

    # Per-cluster average layer assignments.
    cluster_layer_dat <- in_celltype_meta %>%
        group_by(clinical_dx, region) %>%
        slice_sample(n = 20000) %>%
        filter(!cluster_celltype_layer %in% c("pericyte", "ependymal", "t_cell")) %>%
        ungroup %>%
        mutate(clinical_dx = fct_relevel(clinical_dx, "Control")) %>%
        group_by(library_id, clinical_dx, region, cluster_celltype_layer) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(library_id) %>%
        mutate(percent_of_library = n / sum(n)) %>%
        arrange(region, clinical_dx, cluster_celltype_layer)

    write_csv(layer_dat, file.path(out_dir, "celltype_layer_proportions.csv"))
    write_csv(cluster_layer_dat, file.path(out_dir, "cluster_celltype_layer_proportions.csv"))

    pdf(file.path(out_dir, "cluster_celltype_layer_proprtion_dx_test.pdf"), width = 3 * length(unique(cluster_layer_dat$cluster_celltype_layer)), height = 2 * length(unique(cluster_layer_dat$region)))

    my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
    plt <- ggplot(cluster_layer_dat, aes(x = clinical_dx, y = percent_of_library, fill = clinical_dx)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid("region ~ cluster_celltype_layer") +
        geom_jitter(color = "black", alpha = 0.5, position = position_jitter(seed = 1)) +
        stat_compare_means(
            aes(group = clinical_dx, label = round(as_label(..p.format..), 2)),
            comparisons = my_comparisons,
            method = "wilcox.test",
            p.adjust.method = "none", 
            label.y.npc = 0
        ) +
        coord_flip() +
        ggtitle(str_glue("{date()} {commit_id} cell type percentage of library + significance testing between dx groupings"))
    print(plt)
    graphics.off()
}

if (!interactive()) {
    main()
}
