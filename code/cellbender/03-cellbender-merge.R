set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "reticulate")
source("complexheatmap.R")
source("palette.R")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf)

in_cellbender_tb <- "../../analysis/seurat_lchen/cellbender/cellranger_h5.csv"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
SEQ_META <- "/u/project/geschwind/drewse/g_singlecell/nucseq-nd/data/library_ids_plus-md.csv"
CELLBENDER_OUT_DIR <- "../../analysis/seurat_lchen/cellbender/runs"
OUT_DIR <- "../../analysis/seurat_lchen/cellbender/merge/"
    
sc <- import("scanpy", convert = F)
pywarn <- import("warnings", convert = F)
ad <- import("anndata", convert = F)

plot_ts <- function() {str_glue("{date()} {tools::md5sum('03-cellbender-merge.R')}")}

main <- function() {
    dir.create(OUT_DIR)
    cellbender_tb <- read_csv(in_cellbender_tb)
    seq_meta <- read_csv(SEQ_META)
    excludes <- read_xlsx(clusters_exclude_file)

    CELLBENDER_OUT_DIR <- normalizePath(CELLBENDER_OUT_DIR)
    cellbender_tb <- cellbender_tb %>%
        mutate(
            counts_h5 = str_glue("{dirname(molecule_h5)}/raw_feature_bc_matrix.h5"),
            cellbender_out_dir = str_glue("{CELLBENDER_OUT_DIR}/{library_id}"),
            out_h5 = str_glue("{cellbender_out_dir}/{library_id}_cellbender_filtered.h5"),
            out_cells = str_glue("{cellbender_out_dir}/{library_id}_cellbender_cell_barcodes.csv")
        )
    # Fixup cellranger library ids to match final metadata library ids.
    cellbender_f <- cellbender_tb %>%
        filter(file.exists(out_h5)) %>%
        left_join(seq_meta, by = c("library_id" = "agg_mat$library_id")) %>%
        mutate(
            library_id = fix_sequencing_sample_ids(library_id, seq_meta),
            region = fct_recode(region, "preCG" = "PreCG")
        )

    # Merge per region.
    cellbender_region <- cellbender_f %>%
        group_by(region) %>%
        group_nest() %>%
        mutate( 
            merge_out = str_glue("{OUT_DIR}/{region}.rds")
        )

    pwalk(cellbender_region, function(...) {
        cr <- list(...)
        print(cr$merge_out)
        gc()

        sc_list <- map(cr$data$out_h5, read_scanpy_h5)
        names(sc_list) <- cr$data$library_id
        ad_concat <- ad$concat(dict(sc_list), label = "library_id", index_unique = "-")

        ad_mat <- py_to_r(ad_concat$X$transpose())
        colnames(ad_mat) <- py_to_r(ad_concat$obs_names$tolist())
        rownames(ad_mat) <- py_to_r(ad_concat$var_names$tolist())

        merge_sobj <- CreateSeuratObject(ad_mat, meta.data = py_to_r(ad_concat$obs))
        merge_sobj[["UMI"]] <- str_extract(colnames(ad_mat), "\\w+-1")

        saveRDS(merge_sobj, cr$merge_out, compress = F)
    })
    saveRDS(cellbender_region, file.path(OUT_DIR, "cellbender_region.rds"), compress = F)

}

fix_sequencing_sample_ids <- function(sample_id, seq_meta) {
    rename_vec <- setNames(seq_meta$library_id, nm = seq_meta$`agg_mat$library_id`)
    return(as.character(rename_vec[sample_id]))
}

read_scanpy_h5 <- function(path) { 
    writeLines(str_glue("read_scanpy_h5: {path}"))
    with(pywarn$catch_warnings(), {
        pywarn$simplefilter("ignore")
        sc_py <- sc$read_10x_h5(path)
        sc_py$var_names_make_unique()
    })
    return(sc_py)
}

if (!interactive()) main()
