# plot external reference fdr + ora tests using hclust ordering.
liblist <- c("Seurat", "circlize", "scales", "ComplexHeatmap", "patchwork", "tidyverse")
l <- lapply(liblist, function(x) {suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = TRUE))})

FDR_FILE <- normalizePath("../../resources/liger_subcluster_hclust_refmat/fdr_long.csv", mustWork = T)
ORA_FILE <- normalizePath("../../resources/liger_subcluster_hclust_refmat/ora_long.csv", mustWork = T)
HCLUST_F <- "../../analysis/seurat_lchen/liger_subcluster_hier/hclust/dge_enrichment_var_beta_hclust.rds"
OUT_DIR <- "../../analysis/seurat_lchen/liger_subcluster_hier/heatmap_refmat/"

plot_ts <- function() str_glue(tools::md5sum("liger_subcluster_ref.R"), " liger_subcluster_ref.R ", date())


main <- function() {
    fdr_tb <- read_csv(FDR_FILE)
    ora_tb <- read_csv(ORA_FILE)
    hclust_tb <- readRDS(HCLUST_F)

    dir.create(OUT_DIR)

    pwalk(hclust_tb, function(...) {
        cr <- list(...)
        clusters <- cr$hclust$labels
        fdr_f_tb <- fdr_tb %>%
            filter(ct_subcluster %in% clusters) %>%
            complete(ct_subcluster = clusters) %>%
            mutate(fdr_log = log10(fdr))
        ora_f_tb <- ora_tb %>%
            filter(ct_subcluster %in% clusters) %>%
            complete(ct_subcluster = clusters) %>%
            mutate(ora_log = log1p(ora))
        if (!all(is.na(fdr_f_tb$fdr)) & !all(is.na(ora_f_tb$ora))) {
            print(cr$cluster_cell_type)
            out_pdf <- str_glue("{OUT_DIR}/{cr$cluster_cell_type}_ref.pdf")
            fdr_plot <- mk_complexheatmap(fdr_f_tb, cr$hclust, seq_pval_colormap, "ct_subcluster", "fdr_log", "ext_module") %>%
                wrap_heatmap() + 
                plot_annotation(title = str_glue("{cr$cluster_cell_type} log10 fdr"), subtitle = plot_ts())
            ora_plot <- mk_complexheatmap(ora_f_tb, cr$hclust, seq_log_colormap, "ct_subcluster", "ora_log", "ext_module") %>%
                wrap_heatmap() + 
                plot_annotation(title = str_glue("{cr$cluster_cell_type} log1p ora"), subtitle = plot_ts())
                print
            pdf(out_pdf, width = 10, height = 10)
            plot(fdr_plot)
            plot(ora_plot)
            dev.off()
        }
    })
    graphics.off()
}

seq_log_colormap <- function(mat, n = 10) {
    colorRamp2(
        breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = n),
        colors = viridis_pal()(n)
    ) 
}

seq_pval_colormap <- function(mat, n = 10) {
    colorRamp2(
        breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = n),
        colors = rev(viridis_pal()(n))
    ) 
}

div_colormap <- function(mat) {
    colorRamp2(
        breaks = c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )
}

mk_complexheatmap <- function(plot_tb, hclust_cols, colorfunc,
    rows, value, cols) {
    plot_matrix <- pivot_matrix(plot_tb, rows, value, cols)
    plot_matrix <- plot_matrix[, hclust_cols$labels]
    colormap <- colorfunc(plot_matrix)

    heat <- Heatmap(
        plot_matrix,
        name = " ",
        col = colormap,
        na_col = "grey75",
        cluster_columns = hclust_cols,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        height = unit(10, "cm")
    )
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from, fill_na = NA) {
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = all_of(values_from), values_fill = fill_na) %>%
        ungroup

    tb_matrix <- select(tb_pivot, -all_of(rows_from)) %>% as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(all_of(rows_from)), collapse = "|"), .groups = "drop") %>%
        pluck("join_rowname")

    return(tb_matrix)
}

tmp_pdf <- function(w, h) { pdf(tempfile(), w, h) }

wrap_heatmap <- function(ch) {
    heat_gtree <- grid.grabExpr(
        draw(ch),
        device = tmp_pdf,
        width = unit(10, "cm"),
        height = unit(10, "cm")
    )
    return(wrap_plots(heat_gtree))
}


