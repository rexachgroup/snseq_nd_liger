set.seed(0)
liblist <- c("Seurat", "tidyverse")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)

OUT_DIR <- "../../analysis/seurat_lchen/liger_subcluster_lme/dge_export"
SEURAT <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_wk.rds"
GENESETS <- list(
    "AD" = list("igap", "kunkle"),
    "bvFTD" = "ferrari",
    "PSP-S" = "hoglinger"
)

main <- function() {
    dir.create(OUT_DIR)
    subcluster_tbs <- readRDS(SEURAT)
    subcluster_tbs <- subcluster_tbs %>%
        filter(map_lgl(broom_join, is_tibble))
    subcluster_tbs <- mutate(subcluster_tbs, pivot_data = map(broom_join, pivot_dge_data))
    subcluster_tbs <- mutate(subcluster_tbs, magma_up = map(pivot_data, magma_filter_up))
    subcluster_tbs <- mutate(subcluster_tbs, magma_down = map(pivot_data, magma_filter_down))
    cluster_tbs <- subcluster_tbs %>%
        group_by(region, cluster_cell_type) %>%
        group_nest() %>%
        mutate(data = map(data, function(data) { 
            up_bind <- bind_rows(data$magma_up)
            dn_bind <- bind_rows(data$magma_down)
            bind_rows(list(up = up_bind, down = dn_bind), .id = "type") %>%
                group_by(type, dx)
        }))
        
    cluster_tbs <- cluster_tbs %>%
        unnest(data) %>%
        group_by(region, cluster_cell_type, dx, type) %>%
        group_nest

    cluster_tbs <- cluster_tbs %>%
        mutate(
            name = str_glue("liger_subcluster-{region}-{cluster_cell_type}-{dx}-{type}"),
            csv_filepath = file.path(normalizePath(OUT_DIR), str_glue("{name}.csv")),
            spec_filepath = str_glue("../input/annotations/{name}.csv"),
            genesets = GENESETS[dx]
        )
    pwalk(cluster_tbs, function(...) {
        cr <- list(...)
        write_csv(cr$data, cr$csv_filepath)
    })

    cluster_tbs %>%
        select(name = name, annotations = csv_filepath, genesets = genesets) %>%
        saveRDS(file.path(OUT_DIR, "annotation_spec.rds"))
}

magma_filter_up <- function(data) {
    data %>%
        group_by(gene, clinical_dx) %>%
        pivot_wider(names_from = "type", values_from = "value") %>%
        filter(estimate > 0, p.value.adj < 0.05) %>%
        select(GENE = gene, GROUP = ct_subcluster, dx = clinical_dx)
}

magma_filter_down <- function(data) {
    data %>%
        group_by(gene, clinical_dx) %>%
        pivot_wider(names_from = "type", values_from = "value") %>%
        filter(estimate < 0, p.value.adj < 0.05) %>%
        select(GENE = gene, GROUP = ct_subcluster, dx = clinical_dx)
}

pivot_dge_data <- function(data) {
    pivot_spec <- build_longer_spec(
        data,
        cols = starts_with("clinical_dx"),
        names_to = c("clinical_dx", "type"),
        names_pattern = "clinical_dx([^.]+)\\.(.+)",
        values_to = "value"
    )
    data_longer <- pivot_longer_spec(data, pivot_spec)
    return(data_longer)
}

