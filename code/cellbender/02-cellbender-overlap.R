# Compare dropped barcodes against liger subclusters.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "batchtools", "hdf5r", "ComplexHeatmap", "patchwork", "circlize")
source("complexheatmap.R")
source("palette.R")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_seurat_meta <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_cellbender_tb <- "../../analysis/seurat_lchen/cellbender/cellranger_h5.csv"
batchtools <- "../../analysis/seurat_lchen/cellbender/batchtools"
cellbender_env <- "../../ext/cellbender/conda/"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
SEQ_META <- "/u/project/geschwind/drewse/g_singlecell/nucseq-nd/data/library_ids_plus-md.csv"
CELLBENDER_OUT_DIR <- "../../analysis/seurat_lchen/cellbender/runs"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/"

plot_ts <- function() {str_glue("{date()} {tools::md5sum('02-cellbender-overlap.R')}")}

main <- function() {
    cellbender_tb <- read_csv(in_cellbender_tb)
    meta <- readRDS(in_seurat_meta)
    seq_meta <- read_csv(SEQ_META)
    excludes <- read_xlsx(clusters_exclude_file)

    CELLBENDER_OUT_DIR <- normalizePath(CELLBENDER_OUT_DIR)
    cellbender_tb <- cellbender_tb %>%
        mutate(
            counts_h5 = str_glue("{dirname(molecule_h5)}/raw_feature_bc_matrix.h5"),
            out_dir = str_glue("{CELLBENDER_OUT_DIR}/{library_id}"),
            out_h5 = str_glue("{out_dir}/{library_id}_cellbender.h5"),
            out_cells = str_glue("{out_dir}/{library_id}_cellbender_cell_barcodes.csv")
        )
    cellbender_f <- cellbender_tb %>%
        filter(file.exists(out_h5))

    # Read cellbender and get cell probability per barcode.
    cellbender_f <- cellbender_f %>%
        mutate(cellbender_prob_df = map(out_h5, read_cellbender_cell_probability))

    # Select library id, barcode, probability columns.
    # Fixup cellranger library ids to match final metadata library ids.
    cellbender_prob_tb <- cellbender_f %>%
        select(library_id, cellbender_prob_df) %>%
        unnest(cellbender_prob_df) %>%
        mutate(
            barcode = str_extract(barcode, "^\\w+"),
            library_id = fix_sequencing_sample_ids(library_id, seq_meta)
        )

    meta_join <- meta %>%
        mutate(barcode = str_extract(UMI, "^\\w+")) %>%
        left_join(cellbender_prob_tb, by = c("library_id", "barcode"))

    # validate matching by printing NAs after join.
    meta_join %>%
        group_by(library_id) %>%
        summarize(n = n(), na = sum(is.na(latent_cell_probability))) %>%
        print(n = Inf)
 
    # plot distribution per cell type.
    gg_ct <- ggplot(meta_join, aes(x = cell_type, y = latent_cell_probability, fill = cell_type)) +
        geom_point(size = 0, alpha = 1, position = "jitter")
    ggsave(file.path(OUT_DIR, "celltype_violin.jpg"), gg_ct, width = 14, height = 7, dpi = 200, quality = 85)

    # plot fraction of non-ambient cells per subcluster + library id.
    frac_summary <- meta_join %>%
        group_by(library_id, clinical_dx, region, ct_subcluster) %>%
        summarize(
            n = n(), 
            lp = sum(latent_cell_probability > 0.5), 
            nonambient_frac = lp / n, 
            median_prob = median(latent_cell_probability),
            .groups = "drop")

    frac_summary %>%
        group_by(region) %>%
        group_nest(keep = T) %>%
        pwalk(., function(...) {
            cr <- list(...)

            plot_meta <- cr$data %>%
                arrange(clinical_dx, library_id, ct_subcluster)

            out_path <- str_glue("{OUT_DIR}/cellbender_frac_cells_{cr$region}.pdf")
            plot_mat <- pivot_matrix(plot_meta, "library_id", "nonambient_frac", "ct_subcluster")
            is_excluded_subcluster <- rownames(plot_mat) %in% excludes$ct_subcluster
            subcluster_order <- order(is_excluded_subcluster)
            is_excluded_subcluster <- is_excluded_subcluster[subcluster_order]
            plot_mat <- plot_mat[subcluster_order,]
            dx_lvls <- unique(plot_meta$clinical_dx) %>% str_sort %>% fct_relevel("Control") %>% levels
            dx_map <- plot_meta[match(colnames(plot_mat), plot_meta$library_id),]$clinical_dx

            excl_anno <- HeatmapAnnotation(
                is_excluded = is_excluded_subcluster,
                col = list(is_excluded = c("FALSE" = "gray60", "TRUE" = "red")),
                which = "row"
            )
            dx_anno <- HeatmapAnnotation(
                dx = dx_map,
                col = list(dx = setNames(pal_stallion[1:4], dx_lvls)),
                which = "col"
            )

            colormap <- colorRamp2(
                breaks = seq(0, 1, 0.1),
                colors = scales::viridis_pal()(11)
            )

            ch <- Heatmap(
                plot_mat,
                name = "frac_cells",
                col = colormap,
                left_annotation = excl_anno,
                bottom_annotation = dx_anno,
                cluster_rows = F,
                cluster_columns = F,
                row_title = "ct_subcluster",
                column_title = "library_id"
            )
            plot_size <- ch_layout_size(ch)
            plot_grob <- ch_to_grob(ch)
            plot_w <- wrap_plots(plot_grob) + plot_annotation(
                title = str_glue("{cr$region} percentage passed cellbender"),
                subtitle = plot_ts()
            )
            ggsave(out_path, plot_w, width = plot_size[[1]] +  3, height = plot_size[[2]] + 8)
        })
}

# Extract cell probability from cellbender output.
read_cellbender_cell_probability <- function(cellbender_h5_path) {
    h5_file <- H5File$new(cellbender_h5_path, mode = "r")
    # https://github.com/broadinstitute/CellBender/issues/152 for h5 file structure and python loading.
    # latent_cell_probablility is thresholded at > 0.5 to produce cellbender barcode file
    # and filtered version of h5 file.

    cellbender_cell_prob_vec <- h5_file[['matrix/latent_cell_probability']][]
    barcode_kept_indices <- h5_file[["matrix/barcode_indices_for_latents"]][]
    barcodes <- h5_file[["matrix/barcodes"]][]
    
    h5_file$close_all()

    filtered_barcodes <- barcodes[barcode_kept_indices + 1] # R is 1-indexed, python / hdf5 is 0-indexed
    return(data.frame(barcode = filtered_barcodes, latent_cell_probability = cellbender_cell_prob_vec))
}

fix_sequencing_sample_ids <- function(sample_id, seq_meta) {
    rename_vec <- setNames(seq_meta$library_id, nm = seq_meta$`agg_mat$library_id`)
    return(as.character(rename_vec[sample_id]))
}

if (!interactive()) main()
