# Plot per-dx lme betas from liger_subcluster_lme on a heatmap.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "scales", "readxl", "batchtools", "ComplexHeatmap", "circlize")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
options(future.globals.maxSize = Inf, deparse.max.lines = 5)

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
            (region == "preCG" & cluster_cell_type == "opc")
        )
    
    dx <- c("AD", "bvFTD", "PSP-S")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_pval_cols <- str_glue("clinical_dx{dx}.p.value")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    clinical_dx_z_cols <- str_glue("clinical_dx{dx}.statistic")
    filter_args <- tibble(dx, clinical_dx_beta_cols, clinical_dx_pval_cols, clinical_dx_fdr_cols, clinical_dx_z_cols)

    filtered_gene_tests <- full_join(cluster_wk, filter_args, by = character())
    glimpse(filtered_gene_tests)

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
    
    plot_gene_list <- filtered_gene_list %>%
        select(region, dx, cluster_cell_type, ct_subcluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        inner_join(marker_gene_tb, by = "gene") %>%
        mutate(
            gene = fct_rev(fct_relevel(gene, marker_gene_tb$gene)),
            fdr_2 = p.adjust(pval, n = nrow(marker_gene_tb), method = "BH"),
            stars = unclass(symnum(pval, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.005, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " "))),
            ct_subcluster = fct_drop(ct_subcluster),
            column_id = paste(dx, ct_subcluster, sep = "_")
        ) %>%
        glimpse


    #     plot_gene_list_nafill <- plot_gene_list %>%
    #         complete(gene, nesting(region, dx, ct_subcluster))

    plot_matrix_fill <- pivot_matrix(plot_gene_list, "column_id", "z", "gene")
    plot_matrix_text <- pivot_matrix(plot_gene_list, "column_id", "stars", "gene")

    row_annot <- marker_gene_tb %>% filter(gene %in% rownames(plot_matrix_fill))
    plot_row_order <- as.character(row_annot$gene)
    plot_row_split <- as.character(row_annot$disease)
    col_annot <- plot_gene_list %>%
        group_by(dx, ct_subcluster) %>%
        slice_head() %>%
        ungroup
    plot_col_order <- col_annot$column_id
    plot_col_label <- as.character(col_annot$ct_subcluster)
    plot_col_split <- col_annot %>% select(dx, region, cluster_cell_type) #paste(col_annot$dx, col_annot$region, col_annot$cluster_cell_type, sep = "_")

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

    xwidth = 0.25 * ncol(plot_matrix_fill) 
    yheight = 0.25 * nrow(plot_matrix_fill)
    pdf(file.path(out_path_base, "subcluster_lme_heatmap.pdf"), width = xwidth, height = yheight)
    print(Heatmap(
        plot_matrix_fill, 
        col = colormap,
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
    ))
    #     ggplot(plot_gene_list, aes(x = ct_subcluster, y = gene, fill = z, label = stars)) +
    #         geom_tile() +
    #         geom_text() +
    #         facet_wrap(c("dx", "region"), scales = "free", ncol = 99) +
    #         scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), na.value = "grey75") +
    #         theme(axis.text.x = element_text(angle = 90))
    graphics.off()
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

p.adjust <- function (p, method = p.adjust.methods, n = length(p)) {
    method <- match.arg(method)
    if (method == "fdr")
        method <- "BH"
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p)))
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    #stopifnot(n >= lp)
    if (n <= 1)
        return(p0)
    if (n == 2 && method == "hommel")
        method <- "hochberg"
    p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
        i <- seq_len(lp)
        o <- order(p)
        ro <- order(o)
        pmin(1, cummax((n + 1L - i) * p[o]))[ro]
    }, hommel = {
        if (n > lp) p <- c(p, rep.int(1, n - lp))
        i <- seq_len(n)
        o <- order(p)
        p <- p[o]
        ro <- order(o)
        q <- pa <- rep.int(min(n * p/i), n)
        for (j in (n - 1L):2L) {
            ij <- seq_len(n - j + 1L)
            i2 <- (n - j + 2L):n
            q1 <- min(j * p[i2]/(2L:j))
            q[ij] <- pmin(j * p[ij], q1)
            q[i2] <- q[n - j + 1L]
            pa <- pmax(pa, q)
        }
        pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
    }, hochberg = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin((n + 1L - i) * p[o]))[ro]
    }, BH = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }, BY = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        q <- sum(1/(1L:n))
        pmin(1, cummin(q * n/i * p[o]))[ro]
    }, none = p)
    p0
}

if (!interactive()) {
    main()
}
