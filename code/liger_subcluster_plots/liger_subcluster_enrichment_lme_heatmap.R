set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "scales", "lmerTest", "broom.mixed", "future.apply")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_plots/"

in_cluster_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
marker_gene_tb <- tibble(gene = c("KCNH7", "NLGN1", "RORA", "RORB", "PDE1C", "OPCML"))

main <- function() {
    cluster_wk_in <- readRDS(in_cluster_wk)
    excludes <- read_xlsx(clusters_exclude_file)

    cluster_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        filter(
            (region == "insula" & cluster_cell_type == "excitatory") |
            (region == "preCG" & cluster_cell_type == "excitatory") |
            (region == "calcarine" & cluster_cell_type == "excitatory")
        )

    # re-run lmer for just markers.

    model_design <- "expression ~ enrichment_cluster + pmi + age + sex + number_umi + percent_mito + (1 | library_id)"
    
    cluster_wk$marker_expr <- pmap(cluster_wk, function(...) {
            cr <- list(...)
            get_gene_expr(cr$data, cr$ct_subcluster, marker_gene_tb$gene)
        })

    cluster_wk$lm_test <- pmap(cluster_wk, function(...) {
            cr <- list(...)
            marker_expr <- cr$marker_expr
            test_cluster <- cr$liger_clusters
            test_vars <- cr$data %>%
                mutate(
                    enrichment_cluster = fct_other(liger_clusters, keep = test_cluster, other_level = "all_other_cells") %>% fct_relevel("all_other_cells")
                )
            model <- as.formula(model_design)
            run_lmer(marker_expr, test_vars, model, 8)
        })


    cluster_wk$lm_join <- pmap(cluster_wk, function(...) {
            cr <- list(...)
            cr$lm_test %>%
                bind_rows(.id = "gene") %>%
                filter(!is.na(term) & grepl("enrichment_cluster", term)) %>%
                mutate(fdr = p.adjust(p.value, method = "BH"))
        })
    
    # unnest, format, and set up heatmap.
    heatmap_lme_tb <- cluster_wk %>% select(-data, -lm_test, -broom_join, -marker_expr) %>% unnest(lm_join)
    heatmap_lme_tb <- heatmap_lme_tb %>%
        mutate(
            gene = fct_rev(fct_relevel(gene, marker_gene_tb$gene)),
            fdr_2 = p.adjust(p.value, n = nrow(marker_gene_tb), method = "BH"),
            fdr_display = signif(fdr_2, 2),
            stars = unclass(symnum(fdr_2, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.005, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " "))),
            ct_subcluster = fct_drop(ct_subcluster)
        )

    xwidth <- 1 * length(unique(heatmap_lme_tb$ct_subcluster)) + 1
    yheight <- 1 * length(unique(heatmap_lme_tb$gene)) + 2
    pdf(file.path(out_path_base, "subcluster_lme_enrichment_heatmap.pdf"), width = xwidth, height = yheight)
    print(heatmap_text_matrix(heatmap_lme_tb, "ct_subcluster", "gene", "estimate", "fdr_display"))
    graphics.off()

    write_csv(heatmap_lme_tb, file.path(out_path_base, "subcluster_lme_enrichment_tb.csv"))

}

get_gene_expr <- function(meta, ct_subcluster, genes) {
    print(dim(meta))
    envdata <- new.env()
    load(file.path("..", unique(meta$filepath)), envdata)
    print(dim(envdata$liger_so))

    expr_m <- GetAssayData(envdata$liger_so, slot = "data")
    expr_m <- expr_m[genes, meta$cell_ids]
    print(dim(expr_m))
    rm(envdata)
    return(expr_m)
}

run_lmer <- function(expr_m, test_vars, model, cores = NULL) {
    if (!is.null(cores)) {
        plan(multicore, workers = cores)
    }
    lm_out_obj_l <- future_apply(X = expr_m, MARGIN = 1, FUN = function(expr_r) {
        dat <- data.frame(expression = as.vector(expr_r), test_vars)
        # use tryCatch to return NA when model can't be fit for a gene
        tryCatch({
            broom::tidy(lmerTest::lmer(model, data = dat))
        },
        error = function(x) { print(x); return(x) })
    })
    return(lm_out_obj_l)
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

heatmap_text_matrix <- function(tb, row_str, col_str, val_str, txt_str, row_split = NULL, col_split = NULL) {
    heatmap_val_matrix <- pivot_matrix(tb, row_str, val_str, col_str)
    heatmap_text_matrix <- pivot_matrix(tb, row_str, txt_str, col_str)
    row_annot <- tb %>%
        group_by(.data[[col_str]]) %>%
        slice_head %>% ungroup
    col_annot <- tb %>%
        group_by(.data[[row_str]]) %>%
        slice_head %>% ungroup

    plot_row_order <- str_sort(row_annot[[col_str]])
    plot_col_order <- str_sort(col_annot[[row_str]])
    
    plot_row_split <- NULL
    plot_col_split <- NULL
    
    plot_row_label <- as.character(plot_row_order)
    plot_col_label <- as.character(plot_col_order)

    heatmap_val_matrix <- heatmap_val_matrix[plot_row_order, plot_col_order]
    heatmap_text_matrix <- heatmap_text_matrix[plot_row_order, plot_col_order]
    
    plot_text_label <- function(j, i, x, y, w, h, col) {
        label <- heatmap_text_matrix[i, j]
        if (!is.na(label)) {
            grid.text(label, x, y)
        }
    }
    colormap <- colorRamp2(
        breaks = c(min(heatmap_val_matrix, na.rm = TRUE), 0, max(heatmap_val_matrix, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )

    Heatmap(
        heatmap_val_matrix,
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
    )
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
