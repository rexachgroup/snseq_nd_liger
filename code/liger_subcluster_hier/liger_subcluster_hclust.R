# Hierarchically clustering Liger clusters based on top genes.
# Defining top genes as top 100 beta variance genes in enrichment dge

set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "scales")
l <- lapply(liblist, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE)))
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_hier/hclust"

in_cluster_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
excitatory_markers <- read_csv("../../resources/excitatory_markers_20191023.csv")

RESOURCES <- list(
    ncpus = 4,
    memory = 32,
    walltime = 18000,
    partition = "bigmem"
)

main <- function() { 
    if (!dir.exists(out_path_base)) dir.create(out_path_base, recursive = TRUE)

    cluster_wk_in <- readRDS(in_cluster_wk)
    excludes <- read_xlsx(clusters_exclude_file)
    cluster_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster)
        #filter(cluster_cell_type == "excitatory", region == "insula" | region == "preCG")
    #     marker_tb <- excitatory_markers %>%
    #         filter(!is.na(gene_symbol)) %>%
    #         filter(!duplicated(gene_symbol)) %>%
    #         arrange(marker_for, gene_symbol)


    cluster_ct_group <- cluster_wk %>%
        group_by(cluster_cell_type) %>%
        group_nest

    cluster_ct_group <- mutate(cluster_ct_group, 
        beta_hclust_data = map(data, ~tryCatch(cluster_worker(.), error = function(x) {print(x); return(NA)}) )) %>%
        filter(!is.na(beta_hclust_data))

    pwalk(cluster_ct_group, function(...) {
        cr <- list(...)
        hclust_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_var_beta_hclust.pdf"))
        heatmap_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_var_beta_heatmap.pdf"))

        pdf(hclust_path, width = 10, height = 10)
        tryCatch(plot(cr$beta_hclust_dat$hclust, main = "hclust on top 100 variable genes by module signature beta"))
        graphics.off()
        print(hclust_path)

        pdf(heatmap_path, width = 25, height = 10)
        tryCatch(print(cr$beta_hclust_data$heatmap))
        graphics.off()
        print(heatmap_path)

    })

    pwalk(cluster_ct_group, function(...) {
        cr <- list(...)
        rds_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_hclust.rds"))
        saveRDS(list(top_var_genes = cr$beta_hclust_data$top_var_genes, hclust = cr$beta_hclust_data$hclust), rds_path)
        print(rds_path)
    })

    print()

#     # Pull top variable genes.
#     # Plot hclust structure.
#     pdf(file.path(out_path_base, "dge_enrichment_var_beta_hclust.pdf"), width = 10, height = 10)
#     tryCatch(plot(beta_hclust, main = "hclust on top 100 variable genes by module signature beta"))
#     graphics.off()
#     
#     # Plot hclust plus beta expression.
#     pdf(file.path(out_path_base, "dge_enrichment_var_beta_heatmap.pdf"), width = 25, height = 10)
#     heatmap_text_matrix(dge_plot_tb, "gene", "ct_subcluster", "beta", row_hclust = TRUE)
#     graphics.off()
# 
#     pdf(file.path(out_path_base, "dge_enrichment_marker_beta_heatmap.pdf"), width = 10, height = 10)
#     dge_marker_beta_matrix <- enrichment_dge_list %>%
#         filter(gene %in% marker_tb$gene_symbol) %>%
#         arrange(var(beta)) %>%
#         heatmap_text_matrix("gene", "ct_subcluster", "beta", row_hclust = TRUE) %>%
#         print
#     graphics.off()
}

cluster_worker <- function(cluster_wk) { 
    broom_tb_func <- function(...) {
        cr <- list(...)
        beta_col <- str_glue("enrichment_cluster{cr$liger_clusters}.estimate")
        pval_col <- str_glue("enrichment_cluster{cr$liger_clusters}.p.value")
        fdr_col <- str_glue("enrichment_cluster{cr$liger_clusters}.p.value.adj")
        z_col <- str_glue("enrichment_cluster{cr$liger_clusters}.statistic")
        cr$broom_join %>%
            rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
            select(gene, beta, pval, fdr, z)
    }
    enrichment_dge_list <- mutate(cluster_wk, broom_tb = pmap(cluster_wk, broom_tb_func)) %>%
        unnest(broom_tb)

    top_var_genes <- enrichment_dge_list %>%
        group_by(gene) %>%
        filter(length(unique(region)) == length(unique(enrichment_dge_list$region))) %>%
        summarize(beta_var = var(beta)) %>%
        arrange(desc(beta_var)) %>%
        slice_head(n = 100) %>%
        print(n = Inf)
    
    dge_plot_tb <- enrichment_dge_list %>%
        filter(gene %in% top_var_genes$gene) %>%
        arrange(var(beta)) %>%
        group_by(gene)

    dge_beta_matrix <- pivot_matrix(dge_plot_tb, "ct_subcluster", "beta", "gene")
    beta_hclust <- hclust(dist(t(dge_beta_matrix)))
    beta_heatmap <- heatmap_text_matrix(dge_plot_tb, "gene", "ct_subcluster", "beta", row_hclust = beta_hclust, legend_name = "DGE beta")
    return(list(var_genes = top_var_genes, matrix = dge_beta_matrix, hclust = beta_hclust, heatmap = beta_heatmap))
}

pivot_matrix <- function(tb, cols_from, values_from, rows_from) {
    tb_pivot <- tb %>%
        select(all_of(c(cols_from, values_from, rows_from))) %>%
        pivot_wider(names_from = all_of(cols_from), values_from = values_from) %>%
        ungroup

    tb_matrix <- select(tb_pivot, -rows_from) %>% as.matrix()

    rownames(tb_matrix) <- tb_pivot %>% 
        rowwise() %>%
        summarize(join_rowname = paste(c_across(rows_from), collapse = "|"), .groups = "drop") %>%
        pluck("join_rowname")

    glimpse(tb_matrix)
    return(tb_matrix)
}

heatmap_text_matrix <- function(tb, 
        row_str, col_str, val_str, txt_str = NULL, 
        row_split = NULL, col_split = NULL, 
        row_hclust = FALSE, col_hclust = FALSE,
        legend_name = ""
){
    heatmap_val_matrix <- pivot_matrix(tb, row_str, val_str, col_str)
    if (!is.null(txt_str)) {
        heatmap_text_matrix <- pivot_matrix(tb, row_str, txt_str, col_str)
    }
    row_annot <- tb %>%
        group_by(.data[[col_str]]) %>%
        slice_head %>% ungroup
    col_annot <- tb %>%
        group_by(.data[[row_str]]) %>%
        slice_head %>% ungroup
    
    if (is.logical(row_hclust) && !row_hclust) {
        plot_row_order <- row_annot[[col_str]]
        plot_row_label <- as.character(plot_row_order)
    }
    
    if (is.logical(col_hclust) && !col_hclust) {
        plot_col_order <- col_annot[[row_str]]
        plot_col_label <- as.character(plot_col_order)
    }

    
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
        name = legend_name,
        col = colormap,
        na_col = "grey75",
        cluster_rows = row_hclust,
        cluster_columns = col_hclust,
        row_title_rot = 0,
        column_title_rot = 90,
        row_title_gp = gpar(fontsize = 10),
        column_title_gp = gpar(fontsize = 10),
        row_names_side = "left"
    )
}

if (!intearactive()) {
    main()
}
