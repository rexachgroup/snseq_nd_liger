# Complexheatmap showing clustering of liger clusters with:
# 1. violin plot of expression values of selected genes.
# 2. counts of genes / cells per cluster.

# GFAP DPP10
# SLC1A2 SLC2A3 GPC5 SOX5
# HPSE2 CACNB2

set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "RColorBrewer", "scales")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
in_cluster_dge_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
in_seurat_liger <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"
clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
excitatory_markers <- "../../resources/excitatory_layers_20210318.csv"
out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_hier/heatmap"
options(future.globals.maxSize = Inf)

plot_ts <- 

main <- function() {
    dir.create(out_path_base, recursive = TRUE, showWarnings = FALSE)
    seurat_obj <- readRDS(in_seurat_rds)
    cluster_wk_in <- readRDS(in_cluster_dge_wk)
    liger_meta <- readRDS(in_seurat_liger)
    excludes <- read_xlsx(clusters_exclude_file)
    marker_tb <- read_csv(excitatory_markers)
    #annot_gene_list <- unique(marker_tb$gene_symbol)
    #annot_gene_list <- c("GFAP", "DPP10", "SLC1A2", "SLC1A3", "GP5", "SOX5", "HPSE2", "CACNB2")
    annot_gene_list <- c("RORB", "KCNH7", "OPTN", "TMEM106B", "GPC5", "GPC6", "GRM8", "GPC4")
    cluster_dge_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster)
    
    cluster_ct_group <- cluster_dge_wk %>%
        group_by(cluster_cell_type) %>%
        group_nest %>%
        filter(cluster_cell_type == "excitatory")

    cluster_ct_group <- mutate(cluster_ct_group, liger_meta = map(data, ~filter(liger_meta, ct_subcluster %in% .$ct_subcluster)))
    cluster_ct_group <- mutate(cluster_ct_group,
        dge_data = map(data, ~tryCatch(fmt_cluster_dge(.), error = function(x) {print(x); return(NA)}))
    ) %>%
    filter(!is.na(dge_data))

    cluster_ct_group <- mutate(cluster_ct_group, dge_plot_tb = map(dge_data, filter_var_genes, filter_genes_all_regions = FALSE))
    cluster_ct_group <- mutate(cluster_ct_group, violin_data = map(liger_meta, mk_gene_data, seurat_obj, annot_gene_list))
    cluster_ct_group <- mutate(cluster_ct_group, barplot_data = map(liger_meta, mk_barplot_data))
    
    cluster_ct_group <- mutate(cluster_ct_group, dge_base_heatmap = map(dge_plot_tb, mk_beta_heatmap))
    cluster_ct_group <- mutate(cluster_ct_group, violin_annot = map(violin_data, mk_gene_annot))
    cluster_ct_group <- mutate(cluster_ct_group, counts_annot = map(barplot_data, mk_counts_annot))
    cluster_ct_group <- mutate(cluster_ct_group, dge_combo_heatmap = pmap(list(dge_base_heatmap, violin_annot, counts_annot), mk_combo_heatmap))

    pwalk(cluster_ct_group, function(...) {
        cr <- list(...)
        heatmap_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_var_beta_heatmap.pdf"))

        pdf(heatmap_path, width = unit(10, "mm"), height = unit(20, "mm"))
        test <- draw(cr$dge_combo_heatmap, heatmap_legend_list = cr$counts_annot$legends)
        graphics.off()
        print(heatmap_path)
    })
}

fmt_cluster_dge <- function(cluster_wk) {
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
}

filter_var_genes <- function(dge_data, filter_genes_all_regions = FALSE) {
    if (filter_genes_all_regions) {
        gene_dge_list <- dge_data %>%
            group_by(gene) %>% 
            filter(length(unique(region)) == length(unique(dge_data$region)))
    } else {
        gene_dge_list <- dge_data %>%
            group_by(gene) 
    }

    top_var_genes <- gene_dge_list %>%
        summarize(beta_var = var(beta)) %>%
        arrange(desc(beta_var)) %>%
        slice_head(n = 100)
    
    dge_plot_tb <- dge_data %>%
        filter(gene %in% top_var_genes$gene) %>%
        arrange(var(beta)) %>%
        group_by(gene)
    
    return(dge_plot_tb)
}

mk_beta_heatmap <- function(dge_plot_tb, beta_hclust = NULL) { 
    dge_beta_matrix <- pivot_matrix(dge_plot_tb, "ct_subcluster", "beta", "gene")
    if (is.null(beta_hclust)) {
        beta_hclust <- hclust(dist(t(dge_beta_matrix)))
    }
    colormap <- colorRamp2(
        breaks = c(min(dge_beta_matrix, na.rm = TRUE), 0, max(dge_beta_matrix, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )
    heat <- Heatmap(
        dge_beta_matrix,
        name = " ",
        col = colormap,
        na_col = "grey75",
        cluster_columns = beta_hclust,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        height = unit(25, "mm")
    )
    return(heat)
}

# - Subset seurat object to cells in dge_plot_tb / genes in annot_gene_list.
# - For each gene, get expression matrix and convert to tibble with cols expr / umi.
#   - join with cell types' metadata and split along ct_subcluster to get a list of vectors containing each subclusters' expression values.
# Returns a list of lists:
# - names(x) is gene names in annot_gene list.
# - names(x[[1]]) is subclusters in liger_meta.
# - x[[1]] is a vector of gene expression values for each subcluster.
mk_gene_data <- function(liger_meta, sobj, annot_gene_list) {
    feature_list <- intersect(annot_gene_list, rownames(seurat_obj))
    #subset_sobj <- subset(sobj, cells = liger_meta$UMI, features = feature_list)
    sobj_mat <- GetAssayData(sobj, slot = "scale.data")
    feature_list <- intersect(annot_gene_list, rownames(sobj_mat))
    subcluster_expr_vecs <- map(feature_list, function(gene_name) {
        gene_expr_tb <- sobj_mat[gene_name,] %>% as.data.frame %>% rename(gene_expr = '.') %>% rownames_to_column("UMI") 
        meta_expr <- left_join(liger_meta, gene_expr_tb, by = "UMI") 
        meta_nest <- meta_expr %>%
            select(ct_subcluster, gene_expr) %>%
            nest(data = -ct_subcluster) %>%
            mutate(data = map(data, ~as.double(pluck(., 1))))
        subcluster_expr_vec_list <- setNames(meta_nest$data, meta_nest$ct_subcluster)
        return(subcluster_expr_vec_list)
    }) %>% setNames(feature_list)
    return(subcluster_expr_vecs)
}

#   - pass to columnAnnotation to create one violin plot per subcluster.
mk_gene_annot <- function(violin_data) {
    imap(violin_data, function(subcluster_expr_vec_list, gene_name) { 
        columnAnnotation(gene = anno_density(subcluster_expr_vec_list, type = "heatmap"), name = gene_name, annotation_label = gene_name, annotation_height = unit(15, "mm"))
    })
}

mk_barplot_data  <- function(liger_meta) {
    cell_count <- liger_meta %>%
        group_by(ct_subcluster) %>%
        summarize(n = n())
    cell_count_list <- setNames(cell_count$n, cell_count$ct_subcluster)

    cell_dx_count <- liger_meta %>%
        group_by(ct_subcluster, clinical_dx) %>%
        summarize(n = n(), .groups = "drop")
    cell_dx_mat <- pivot_matrix(cell_dx_count, "clinical_dx", "n", "ct_subcluster", fill_na = 0)

    gene_count <- liger_meta %>%
        select(ct_subcluster, number_genes) %>%
        nest(data = -ct_subcluster) %>%
        mutate(data = map(data, ~as.double(pluck(., 1))))
    gene_count_list <- setNames(gene_count$data, gene_count$ct_subcluster)

    return(list(cell_count = cell_count_list, cell_dx_count = cell_dx_mat, gene_count = gene_count_list))
}

pal_stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

# make annotations based on cell counts.
# adding legends using annotation_legend_param is not supported -- return as list to plot later down the line
mk_counts_annot <- function(count_data) {
    dx_names <- colnames(count_data$cell_dx_count)
    dx_annot_colors <- pal_stallion[1:length(dx_names)]
    return(
        list(
            annotations = list(
                count_annot = columnAnnotation(
                    count = anno_barplot(count_data$cell_count),
                    annotation_label = "# cells",
                    annotation_height = unit(15, "mm"),
                    annotation_name_rot = 90
                ),
                count_dx_annot = columnAnnotation(
                    count_dx = anno_barplot(count_data$cell_dx_count, gp = gpar(fill = dx_annot_colors)), 
                    annotation_label = "cells\nper\ndx",
                    annotation_height = unit(20, "mm"),
                    annotation_name_rot = 90 
                ),
                gene_annot = columnAnnotation(
                    gene = anno_boxplot(count_data$gene_count, outline = FALSE),
                    annotation_label = "genes\nper\ncell",
                    annotation_height = unit(15, "mm"),
                    annotation_name_rot = 90
                )
            ),
            legends = list(
                count_dx_annot = Legend(labels = dx_names, title = "dx", legend_gp = gpar(fill = dx_annot_colors))
            )
        )
    )
}

mk_combo_heatmap <- function(base_heatmap, violin_annots, counts_annot) {
    hmap <- base_heatmap
    for (annot in violin_annots) {
        hmap <- hmap %v% annot
    }
    for (annot in counts_annot$annotations) {
        hmap <- hmap %v% annot
    }
    # names are dropped after adding annotations. recover by manually adding text annotation.
    # col_labels are reordered during heatmap draw -- don't reorder manually
    base_col_hclust <- base_heatmap@column_dend_param
    base_col_labels <- base_col_hclust$obj$labels
    hmap <- hmap %v% columnAnnotation(colname_annot = anno_text(base_col_labels, location = unit(1, "npc"), just = "right"))
    return(hmap)
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

