# Per-dx celltype proportions for all control samples.
# Downsampling per dx / region to get consistent propoertion estimates(?) between different regions.
liblist <- c("tidyverse", "Seurat", "ggpubr", "broom", "rstatix")
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
        group_by(library_id, autopsy_id, clinical_dx, region, cluster_celltype_layer) %>%
        summarise(n = n(), percent_mito = mean(percent_mito), age = first(age), pmi = first(pmi), sex = first(sex), .groups = "drop") %>%
        group_by(library_id) %>%
        mutate(percent_of_library = n / sum(n)) %>%
        arrange(region, clinical_dx, cluster_celltype_layer)
    
    meta_lm <- cluster_layer_dat %>%
        group_by(cluster_celltype_layer, region) %>%
        summarize(lm_matrix = list(tidy(lm(percent_of_library ~ clinical_dx + sex + age + pmi + percent_mito), data = .))) %>%
        unnest(lm_matrix) %>%
        filter(grepl("clinical_dx", term))

    meta_lm$group1 = str_replace(meta_lm$term, "clinical_dx", "")
    meta_lm$group2 = "Control"
    meta_lm$y.position = "1"

    meta_lm <- cluster_layer_dat %>%
        group_by(cluster_celltype_layer, region) %>%
        tukey_hsd(percent_of_library ~ clinical_dx) %>%
        add_xy_position(x = "clinical_dx", fun = "mean_se")

    write_csv(layer_dat, file.path(out_dir, "celltype_layer_proportions.csv"))
    write_csv(cluster_layer_dat, file.path(out_dir, "cluster_celltype_layer_proportions.csv"))

    pdf(file.path(out_dir, "cluster_celltype_layer_proprtion_dx_test.pdf"), width = 3 * length(unique(cluster_layer_dat$cluster_celltype_layer)), height = 2 * length(unique(cluster_layer_dat$region)))

  #  my_comparisons <- list(c(1, 2), c(1, 3), c(1, 4))
    #     plt1 <- ggplot(cluster_layer_dat, aes_string(x = "clinical_dx", y = "percent_of_library", color = "clinical_dx")) +
    #         facet_grid(c("region", "cluster_celltype_layer"), scales = "free") +
    #         geom_boxplot() +
    #         geom_jitter(color = "black", alpha = 0.5, position = position_jitter(seed = 1)) +
    #         stat_pvalue_manual(
    #             meta_lm,
    #             label = "p.value"
    #         ) +
    #         coord_flip() +
    #         ggtitle(str_glue("{date()} {commit_id} cell type percentage of library + significance testing between dx groupings"))
    plt2 <- ggboxplot(cluster_layer_dat, x = "clinical_dx", y = "percent_of_library", color = "clinical_dx", facet = c("region", "cluster_celltype_layer")) +
        stat_pvalue_manual(
            meta_lm
        ) +
        coord_flip() +
        ggtitle(str_glue("{date()} {commit_id} cell type percentage of library + significance testing between dx groupings"))
    #     print(plt1)
    print(plt2)
    graphics.off()
}

if (!interactive()) {
    main()
}
