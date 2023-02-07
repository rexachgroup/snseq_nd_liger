# Plot per-dx lme betas from liger_subcluster_lme on a heatmap.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "scales", "readxl", "batchtools", "ComplexHeatmap", "circlize", "patchwork")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)  
plot_ts <- function() str_glue("{system('md5sum liger_subcluster_lme_heatmap.R', intern = T)} {date()}")

in_cluster_wk <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_wk.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
dx_marker_genes <- "../../resources/disease_gene_markers_refined.xlsx"
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"

main <- function() {
    cluster_wk_in <- readRDS(in_cluster_wk)
    excludes <- read_xlsx(clusters_exclude_file)
    marker_gene_tb <- read_xlsx(dx_marker_genes) %>% 
        filter(!duplicated(gene)) %>%
        arrange(disease, gene)
    
    cluster_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        filter(
            cluster_cell_type == "microglia" | 
            cluster_cell_type == "oligodendrocyte" |
            (region == "preCG" & cluster_cell_type == "excitatory") |
            (region == "calcarine" & cluster_cell_type == "excitatory") |
            (region == "calcarine" & cluster_cell_type == "endothelia") |
            (region == "preCG" & cluster_cell_type == "endothelia") |
            (region == "preCG" & cluster_cell_type == "astrocyte") |
            (region == "calcarine" & cluster_cell_type == "astrocyte") |
            (region == "preCG" & cluster_cell_type == "opc") |
            (region == "insula" & cluster_cell_type == "excitatory")
        )
    
    dx <- c("AD", "bvFTD", "PSP-S")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_pval_cols <- str_glue("clinical_dx{dx}.p.value")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    clinical_dx_z_cols <- str_glue("clinical_dx{dx}.statistic")
    filter_args <- tibble(dx, clinical_dx_beta_cols, clinical_dx_pval_cols, clinical_dx_fdr_cols, clinical_dx_z_cols)

    filtered_gene_tests <- full_join(cluster_wk, filter_args, by = character())
    glimpse(filtered_gene_tests)
    
    # rename cols from each dx's broom result to beta, pval, fdr, z.
    filter_call <- function(...) {
        cr <- list(...)
        beta_col <- cr$clinical_dx_beta_cols
        pval_col <- cr$clinical_dx_pval_cols
        fdr_col <- cr$clinical_dx_fdr_cols
        z_col <- cr$clinical_dx_z_cols
        writeLines(str_glue("{cr$ct_subcluster}:  {beta_col} {pval_col} {fdr_col} {z_col}"))
        if (!is.na(cr$broom_join) && nrow(cr$broom_join) > 1 && all(c(beta_col, fdr_col, z_col) %in% colnames(cr$broom_join))) {
            broom_filter <- cr$broom_join %>%
                dplyr::filter(gene %in% marker_gene_tb$gene) %>%
                rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
                select(gene, beta, pval, fdr, z)
        }
    }
    filtered_gene_list <- filtered_gene_tests %>%
        mutate(lme_marker_estimates = pmap(., filter_call)) %>%
        glimpse
    
    # Filter down to marker gene list. Add fdr and stars columns.
    plot_gene_list <- filtered_gene_list %>%
        select(region, dx, cluster_cell_type, ct_subcluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        inner_join(marker_gene_tb, by = "gene") %>%
        group_by(region, cluster_cell_type) %>%
        mutate(
            n = length(unique(gene)) * length(unique(ct_subcluster)),
            gene = fct_rev(fct_relevel(gene, marker_gene_tb$gene)),
            p_val_adj = p.fdr(pval, n = n),
            ct_subcluster = fct_drop(ct_subcluster),
            column_id = paste(dx, ct_subcluster, sep = "_")
        )

    plot_gene_nest <- plot_gene_list %>%
        group_nest(keep = T) %>%
        mutate(heatmap = map(data, plot_lme_result, row_annot))

    plot_gene_nest <- plot_gene_nest %>%
        mutate(        
            xwidth = map(data, ~0.25 * length(unique(.x$column_id))),
            yheight = map(data, ~0.25 * length(unique(.x$gene)))
        )

    pwalk(plot_gene_nest, function(...) {
        cr <- list(...)
        pdf(file.path(out_path_base, str_glue("subcluster_lme_heatmap_{cr$region}_{cr$cluster_cell_type}.pdf")), width = cr$xwidth, height = cr$yheight)
        print(cr$heatmap)
        dev.off()
        write.csv(cr$data, file.path(out_path_base, str_glue("subcluster_lme_heatmap_{cr$region}_{cr$cluster_cell_type}.csv")))
    })
    graphics.off()
}

plot_lme_result <- function(lme_result, row_annot) { 
    lme_result <- ungroup(lme_result)
    # Add stars; save legend for printing
    plot_matrix_stars <- symnum(lme_result$p_val_adj, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", " "))
    lme_result$stars <- as.character(plot_matrix_stars)
    plot_matrix_stars_legend <- as.character(attr(plot_matrix_stars, "legend"))

    # Begin ComplexHeatmap formatting
    plot_matrix_fill <- pivot_matrix(lme_result, "column_id", "z", "gene")
    plot_matrix_text <- pivot_matrix(lme_result, "column_id", "stars", "gene")

    row_annot <- marker_gene_tb %>% filter(gene %in% rownames(plot_matrix_fill))
    plot_row_order <- as.character(row_annot$gene)
    plot_row_split <- as.character(row_annot$disease)
    col_annot <- lme_result %>%
        group_by(dx, ct_subcluster) %>%
        slice_head() %>%
        ungroup
    plot_col_order <- col_annot$column_id
    plot_col_label <- as.character(col_annot$ct_subcluster)
    plot_col_split <- col_annot %>% select(dx, region, cluster_cell_type)
    plot_matrix_fill <- plot_matrix_fill[plot_row_order, plot_col_order]
    plot_matrix_text <- plot_matrix_text[plot_row_order, plot_col_order]

    plot_text_label <- function(j, i, x, y, w, h, col) {
        label <- plot_matrix_text[i, j]
        if (!is.na(label)) {
            grid.text(label, x, y)
        }
    }

    colormap <- colorRamp2(
        breaks = c(min(plot_matrix_fill, na.rm = TRUE), 0, max(plot_matrix_fill, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )

    heatmap <- Heatmap(
        plot_matrix_fill, 
        col = colormap,
        name = "subcluster lme beta",
        na_col = "grey75",
        row_order = plot_row_order,
        column_order = plot_col_order,
        row_split = plot_row_split,
        column_split = plot_col_split,
        column_label = plot_col_label,
        cell_fun = plot_text_label,
        row_title_rot = 0,
        column_title_rot = 90,
        row_title_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10),
        row_names_side = "left"
    )
    
    pdf(tempfile())
    hgrob <- grid.grabExpr(draw(heatmap))
    dev.off()
    
    wrap_plots(hgrob) + plot_annotation(title = str_glue("subcluster lme z \n adj_p_val (n = {lme_result$n}) \n (stars = {plot_matrix_stars_legend})"), subtitle = plot_ts())
}

p.fdr <- function(p, n) {
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p)))
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    if (n <= 1)
        return(p0)
    p0[nna] <- {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }
    p0
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(matches(paste0(c(cols_from, values_from, rows_from), collapse = "|"))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = values_from) %>%
        glimpse
    tb_matrix <- select(tb_pivot, -rows_from) %>%
        as.matrix()
    rownames(tb_matrix) <- tb_pivot %>% rowwise() %>% summarize(paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>% pluck(1)
    return(tb_matrix)
}

if (!interactive()) {
    main()
}
