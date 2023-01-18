liblist <- c("tidyverse", "broom", "edgeR", "limma", "readxl", "ggpubr")
lapply(liblist, require, character.only = TRUE)

SOBJ_FILE <- "../../analysis/pci_import/pci_seurat.rds"
META_FILE <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
OUT_DIR <- "../../analysis/seurat_lchen/subcluster_composition/boxplot/"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
SUBCLUSTER_FILTER_FILE <- "../../analysis/seurat_lchen/liger_subcluster_filtered_props.rds"

main <- function() {
    if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    meta_tb <- readRDS(META_FILE)
    excludes <- read_xlsx(clusters_exclude_file)
    percluster_tb <- readRDS(SUBCLUSTER_FILTER_FILE) %>%
        filter(subcluster_3libs_above_10umi) %>%
        mutate(ct_subcluster = paste(region, ct_subcluster, sep = "-"))

    meta <- meta_tb %>%
        mutate(ct_subcluster = paste(region, cluster_cell_type, liger_clusters, sep = "-"),
            clinical_dx = fct_relevel(clinical_dx, "Control"))

    ctl_filter <- meta %>% 
        group_by(ct_subcluster) %>% 
        summarize(ctl = any(clinical_dx == "Control")) %>% 
        filter(ctl == TRUE)

    library_subcluster_counts_raw <- meta %>%
        group_by(region, library_id, ct_subcluster) %>%
        arrange(region, liger_clusters) %>%
        summarize(cell_type = unique(cell_type), 
                  liger_cluster = unique(liger_clusters),
                  subcluster_ct = n())

    # Count the number of cells for each subcluster that has control samples and is not in ct_subcluster.
    library_subcluster_counts <- meta %>%
        filter(cluster_cell_type == cell_type,
               ct_subcluster %in% percluster_tb$ct_subcluster,
               ct_subcluster %in% ctl_filter$ct_subcluster) %>%
        group_by(region, library_id, ct_subcluster) %>%
        arrange(region, liger_clusters) %>%
        summarize(cell_type = unique(cell_type), 
                  liger_cluster = unique(liger_clusters),
                  subcluster_ct = n())

    # Count the total number of cells + summarize metadata per library + Seurat cluster subset.
    library_celltype_counts_full <- meta %>%
        filter(cluster_cell_type == cell_type) %>%
        group_by(region, library_id, cell_type) %>%
        summarize(
            clinical_dx = unique(clinical_dx),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            cluster_ct = n(),
        )

    # Join tables and pivot to library x subcluster cell counts matrix per Seurat cluster.
    library_pct_raw <- inner_join(library_subcluster_counts_raw, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
        mutate(subcluster_pct_norm = subcluster_ct / cluster_ct)

    library_pct <- inner_join(library_subcluster_counts, library_celltype_counts_full, by = c("region", "library_id", "cell_type")) %>%
        mutate(subcluster_pct_norm = subcluster_ct / cluster_ct)

    library_raw_patchwork <- library_pct_raw %>%
        group_by(ct_subcluster) %>%
        group_split() %>%
        map(function(data) {
            ggplot_subcluster_library_counts(data) +
            ggtitle(str_glue("{unique(data$ct_subcluster)} subcluster count, all liger subclusters"))
        })

    library_ct_patchwork <- library_pct %>%
        group_by(ct_subcluster) %>%
        group_split() %>%
        map(function(data) {
            ggplot_subcluster_library_counts(data) +
            ggtitle(str_glue("{unique(data$ct_subcluster)} subcluster count, control + umi filter"))
        })
    
    library_pct_patchwork <- library_pct %>%
        group_by(ct_subcluster) %>%
        group_split() %>%
        map(function(data) {
            ggplot_subcluster_library_counts(data, ct_col = "subcluster_pct_norm") +
            ggtitle(str_glue("{unique(data$ct_subcluster)} subcluster / library pct, control + umi filter"))
        })


    pdf(file.path(OUT_DIR, "counts_raw_boxplot.pdf"), width = 10, height = 10)
    print(library_raw_patchwork)
    graphics.off()

    pdf(file.path(OUT_DIR, "counts_boxplot.pdf"), width = 10, height = 10)
    print(library_ct_patchwork)
    graphics.off()

    pdf(file.path(OUT_DIR, "pct_boxplot.pdf"))
    print(library_pct_patchwork)
    graphics.off()

}

# data$clinical_dx, data$subcluster_count, data$ct_subcluster
ggplot_subcluster_library_counts <- function(data, ct_col = "subcluster_ct") {
    comparisons <- list(c("Control", "AD"), c("Control", "bvFTD"), c("Control", "PSP-S"))
    max_dat <- quantile(data[[ct_col]], 0.95)
    range <- max(data[ct_col]) - min(data[[ct_col]])
    comparison_y <- c(max_dat, max_dat + (0.05 * range), max_dat + (0.10 * range))
    ggplot(data, aes_string(x = "clinical_dx", y = ct_col, color = "clinical_dx")) +
        geom_boxplot() +
        geom_jitter() +
        stat_compare_means(data = data,comparisons = comparisons, label.y = comparison_y)
}

if (!interactive()) { main() }
