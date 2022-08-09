# recreate subcluster_wk from dge csv results.
set.seed(0)
liblist <- c("Seurat", "tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

CSV_DIR <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/enrichment_tables"
OUT_DIR <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/dge_export"

main <- function() {
    dir.create(OUT_DIR)
    csv_tb <- tibble(path = list.files(CSV_DIR, full.names = TRUE))
    csv_tb <- mutate(csv_tb, data = map(path, read_csv))
    csv_tb <- mutate(csv_tb, pivot_data = map(data, pivot_enrichment_data))
    csv_tb <- mutate(csv_tb, magma_data = map(pivot_data, magma_filter_data))
    csv_tb <- mutate(csv_tb, 
        csv_name = str_match(basename(path), "(.*)\\.csv")[,2],
        csv_basename = str_glue("{csv_name}-magma.csv", .sep = ''),
        out_csv = file.path(OUT_DIR, csv_basename)
    )
    pmap(csv_tb, function(...) {
        cr <- list(...)
        write_csv(cr$magma_data, cr$out_csv)
    })
    csv_tb %>%
        select(csv_name, csv_basename) %>%
        mutate(annotations = file.path("../input/genesets", csv_basename)) %>%
        select(name = csv_name, annotations) %>%
        write_csv(file.path(OUT_DIR, "spec.csv"))

    repivot_tb <- repivot_data_cluster_wk(csv_tb)
    repivot_data_ldsc(csv_tb) %>%
        write_tsv(file.path(OUT_DIR, "liger_ldsc.tsv"), col_names = FALSE)
}

pivot_enrichment_data <- function(data) {
    pivot_spec <- build_longer_spec(
        data,
        cols = starts_with("enrichment_cluster"),
        names_to = c("enrichment_cluster", "type"),
        names_pattern = "enrichment_cluster(\\d+)\\.(.+)",
        values_to = "value"
    )
    data_longer <- pivot_longer_spec(data, pivot_spec)
    return(data_longer)
}

magma_filter_data <- function(data) {
    data %>%
        filter(type == "estimate") %>%
        filter(value > 0.3) %>%
        arrange(enrichment_cluster) %>%
        select(GROUP = enrichment_cluster, GENE = gene)
}

repivot_data_cluster_wk <- function(csv_tb) {
    tb <- csv_tb %>%
        select(name = csv_name, data = pivot_data) %>%
        unnest(data)
    tb_nest <- tb %>%
        group_by(region, cell_type, enrichment_cluster) %>%
        group_nest(keep = TRUE)
    tb_nest <- tb_nest %>%
        mutate(ct_subcluster = str_glue("{region}-{cell_type}-{enrichment_cluster}"))  %>%
        mutate(ct_subcluster = fct_relevel(ct_subcluster, str_sort(ct_subcluster, numeric = TRUE)))
    
    tb_nest <- tb_nest %>%
        mutate(cluster_wk = map(data, function(data) {
            pivot_wider(data, names_from = c("enrichment_cluster", "type"), values_from = "value", names_glue = "enrichment_cluster{enrichment_cluster}.{type}")
        }))

    tb_nest <- tb_nest %>%
        rename(liger_clusters = enrichment_cluster)
    return(tb_nest)
}

repivot_data_ldsc <- function(csv_tb) {
    tb <- csv_tb %>%
        select(name = csv_name, data = pivot_data) %>%
        unnest(data)

    tb %>%
        filter(type == "estimate", value > 0.3) %>%
        mutate(ct_subcluster = str_glue("{region}-{cell_type}-{enrichment_cluster}")) %>%
        select(gene = gene, cluster = ct_subcluster)
}
