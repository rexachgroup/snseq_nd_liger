# Plot liger cluster gene enrichment beta from liger_subcluster_enrichment_dge on a heatmap.
set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "scales")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"

in_cluster_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"

main <- function() {
    cluster_wk_in <- readRDS(in_cluster_wk)
    excludes <- read_xlsx(clusters_exclude_file)
    marker_gene_tb <- tibble(gene = c("KCNH7", "NLGN1", "RORA", "RORB", "PDE1C", "OPCML"))

    cluster_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        filter(
            (region == "preCG" & cluster_cell_type == "excitatory") |
            (region == "calcarine" & cluster_cell_type == "excitatory")
        )

    glimpse(cluster_wk)
    
    filter_call <- function(...) {
        cr <- list(...)
        if (!is.na(cr$broom_join) && nrow(cr$broom_join) >= 1) {
            beta_col <- str_subset(colnames(cr$broom_join), "estimate$")
            pval_col <- str_subset(colnames(cr$broom_join), "p.value$")
            fdr_col <- str_subset(colnames(cr$broom_join), "p.value.adj$")
            z_col <- str_subset(colnames(cr$broom_join), "statistic$")
            writeLines(str_glue("{cr$ct_subcluster}:  {beta_col} {pval_col} {fdr_col} {z_col}"))
            broom_filter <- cr$broom_join %>%
                dplyr::filter(gene %in% marker_gene_tb$gene) %>%
                rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
                select(gene, beta, pval, fdr, z)
        }
    }
    filtered_gene_list <- cluster_wk %>%
        mutate(lme_marker_estimates = pmap(., filter_call)) %>%
        glimpse
    
    plot_gene_list <- filtered_gene_list %>%
        select(region, cluster_cell_type, ct_subcluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        left_join(marker_gene_tb, by = "gene") %>%
        mutate(
            gene = as.character(gene, marker_gene_tb$gene),
            fdr_2 = p.adjust(pval, n = nrow(marker_gene_tb), method = "BH"),
            stars = unclass(symnum(fdr_2, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.005, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " "))),
            ct_subcluster = as.character(ct_subcluster)
        ) %>%
        glimpse
    
    plot_matrix_fill <- pivot_matrix(plot_gene_list, "ct_subcluster", "z", "gene")
    plot_matrix_text <- pivot_matrix(plot_gene_list, "ct_subcluster", "stars", "gene")

    row_annot <- plot_gene_list %>%
        group_by(gene) %>%
        slice_head() %>%
        ungroup
    plot_row_order <- row_annot$gene
    plot_row_label <- row_annot$gene
    plot_row_split <- NULL#row_annot %>% select(gene)

    col_annot <- plot_gene_list %>%
        group_by(ct_subcluster) %>%
        slice_head() %>%
        ungroup
    plot_col_order <- as.character(col_annot$ct_subcluster)
    plot_col_label <- as.character(col_annot$ct_subcluster)
    plot_col_split <- as.character(col_annot$region)

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

    xwidth = 0.5 * ncol(plot_matrix_fill) 
    yheight = 1 * nrow(plot_matrix_fill)
    pdf(file.path(out_path_base, "subcluster_lme_gene_heatmap.pdf"), width = xwidth, height = yheight)
    print(Heatmap(
        plot_matrix_fill, 
        col = colormap,
        na_col = "grey75",
        row_order = plot_row_order,
        row_labels = plot_row_label,
        row_split = plot_row_split,
        column_order = plot_col_order,
        column_label = plot_col_label,
        column_split = plot_col_split,
        cell_fun = plot_text_label,
        row_title_rot = 0,
        column_title_rot = 90,
        row_title_gp = gpar(fontsize = 12),
        row_names_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 12),
        row_names_side = "left"
    ))
    graphics.off()
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(matches(paste0(c(cols_from, values_from, rows_from), collapse = "|"))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = values_from)
    tb_matrix <- select(tb_pivot, -rows_from) %>%
        as.matrix()
    print(colnames(tb_matrix))
    rownames(tb_matrix) <- tb_pivot %>% rowwise() %>% summarize(paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>% pluck(1)
    print(rownames(tb_matrix))
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
