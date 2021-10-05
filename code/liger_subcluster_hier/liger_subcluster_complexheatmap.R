# Complexheatmap showing clustering of liger clusters with:
# 1. violin plot of expression values of selected genes.
# 2. counts of genes / cells per cluster.

set.seed(0)
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "RColorBrewer", "scales", "WGCNA")
l <- lapply(liblist, require, character.only = TRUE, quietly = TRUE)
in_cluster_dge_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
in_limma <- "../../analysis/seurat_lchen/subcluster_composition/limma/subcluster_composition_dx_pmi.rds"
in_seurat_liger <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"

clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
excitatory_markers <- "../../resources/excitatory_layers_20210318.csv"
in_bulk_meta <- "../../resources/individual_nd.rds"

out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_hier/heatmap"
batchtools <- file.path(out_path_base, "batchtools")
options(future.globals.maxSize = Inf)

plot_ts <- str_glue("{system('md5sum liger_subcluster_complexheatmap.R', intern = TRUE)}")

main <- function() {
    dir.create(out_path_base, recursive = TRUE, showWarnings = FALSE)

    seurat_obj <- readRDS(in_seurat_rds)
    cluster_wk_in <- readRDS(in_cluster_dge_wk)
    liger_meta <- readRDS(in_seurat_liger)
    bulk_meta <- readRDS(in_bulk_meta)
    excludes <- read_xlsx(clusters_exclude_file)
    limma_tb <- readRDS(in_limma)
    #marker_tb <- read_csv(excitatory_markers)
    #annot_gene_list <- unique(marker_tb$gene_symbol)
    #annot_gene_list <- c("GFAP", "DPP10", "SLC1A2", "SLC1A3", "GP5", "SOX5", "HPSE2", "CACNB2")
    annot_gene_list <- c("RORB", "KCNH7", "OPTN", "TMEM106B", "GPC5", "GPC6", "GRM8", "GPC4")
    bulk_meta <- mutate(bulk_meta, Autopsy.ID = paste0("P", Autopsy.ID))
    nd_tb <- bulk_meta %>%
        group_by(Autopsy.ID, type) %>%
        slice_head(n = 1) %>%
        select(Autopsy.ID, type, score)

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

    cluster_ct_group <- mutate(cluster_ct_group, dge_hclust = map(dge_plot_tb, mk_dge_hclust))
    cluster_ct_group <- mutate(cluster_ct_group, dge_plot_tb = map(dge_data, filter_var_genes, filter_genes_all_regions = FALSE))
    cluster_ct_group <- mutate(cluster_ct_group, gene_expr_data = map(liger_meta, mk_gene_data, seurat_obj, annot_gene_list))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_counts_data = map(liger_meta, mk_cluster_counts_data))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_limma_data = pmap(list(dge_data, dge_hclust), mk_cluster_limma_data, limma_tb))
    cluster_ct_group <- mutate(cluster_ct_group, nd_correlation = map(liger_meta, mk_nd_data, seurat_obj, nd_tb, annot_gene_list))
    
    cluster_ct_group <- mutate(cluster_ct_group, dge_base_heatmap = pmap(list(dge_plot_tb, dge_hclust), mk_beta_heatmap))
    cluster_ct_group <- mutate(cluster_ct_group, gene_expr_annot = map(gene_expr_data, mk_gene_annot))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_counts_annot = map(cluster_counts_data, mk_cluster_counts_annot))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_limma_annot = map(cluster_limma_data, mk_cluster_limma_annot))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_nd_heatmap = map(nd_correlation, mk_nd_heatmap))
    cluster_ct_group <- mutate(cluster_ct_group, dge_combo_heatmap = pmap(cluster_ct_group, mk_combo_heatmap))

    pwalk(cluster_ct_group, function(...) {
        cr <- list(...)
        heatmap_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_var_beta_heatmap.pdf"))

        pdf(heatmap_path, width = unit(15, "mm"), height = unit(20, "mm"))
        test <- draw(cr$dge_combo_heatmap, heatmap_legend_list = c(cr$cluster_counts_annot$legends, cr$cluster_limma_annot$legends))
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

mk_dge_hclust <- function(dge_plot_tb) {
    dge_beta_matrix <- pivot_matrix(dge_plot_tb, "ct_subcluster", "beta", "gene")
    hclust(dist(t(dge_beta_matrix)))
}

mk_beta_heatmap <- function(dge_plot_tb, beta_hclust) { 
    dge_beta_matrix <- pivot_matrix(dge_plot_tb, "ct_subcluster", "beta", "gene")
    dge_beta_matrix[, beta_hclust$labels]
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
# - Regress out variables for display purposes using negative binomial model.
# - For each gene, get expression matrix and convert to tibble with cols expr / umi.
#   - join with cell types' metadata and split along ct_subcluster to get a list of vectors containing each subclusters' expression values.
# Returns a list of lists:
# - names(x) is gene names in annot_gene_list.
# - names(x[[*]]) is subclusters in liger_meta.
# - x[[*]] is a vector of gene expression values for each subcluster.
mk_gene_data <- function(liger_meta, sobj, annot_gene_list) {
    feature_list <- intersect(annot_gene_list, rownames(seurat_obj))
    subset_sobj <- subset(sobj, cells = liger_meta$UMI, features = feature_list)
    #sobj_scaled <- ScaleData(subset_sobj, model.use = "negbinom", vars.to.regress = c("pmi", "age", "sex", "number_umi", "percent_mito"), assay = "RNA")
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

# Run bicor of nd values against gene markers per subcluster.
mk_nd_data <- function(liger_meta, sobj, nd_tb, annot_gene_list) {
    feature_list <- intersect(annot_gene_list, rownames(GetAssayData(sobj, slot = "scale.data")))
    subset_sobj <- subset(sobj, cells = liger_meta$UMI, features = feature_list)
    sobj_mat <- GetAssayData(sobj, slot = "scale.data")

    subcluster_cor_tbs <- map(feature_list, function(gene_name) {
        gene_expr_tb <- sobj_mat[gene_name,] %>% as.data.frame %>% rename(gene_expr = '.') %>% rownames_to_column("UMI") 
        meta_expr <- left_join(liger_meta, gene_expr_tb, by = "UMI") 
        meta_nd_tb <- left_join(meta_expr, nd_tb, by = c("autopsy_id" = "Autopsy.ID"))

        meta_nd_tb %>%
            group_by(ct_subcluster) %>%
            summarize(bicor = bicor(gene_expr, score, use = "pairwise.complete.obs")[, 1]) %>%
            rename({{ gene_name }} := bicor)
    }) 
    return(subcluster_cor_tbs)

}

#   - pass to columnAnnotation to create one violin plot per subcluster.
mk_gene_annot <- function(gene_data) {
    annotation_height <- unit(15, "mm")
    heatmaps <- imap(gene_data, function(subcluster_expr_vec_list, gene_name){
        columnAnnotation(gene = anno_density(subcluster_expr_vec_list, type = "heatmap"), annotation_label = gene_name, annotation_height = annotation_height)
    })
    violins <- imap(gene_data, function(subcluster_expr_vec_list, gene_name){
        columnAnnotation(gene = anno_density(subcluster_expr_vec_list, type = "violin"), annotation_label = gene_name, annotation_height = annotation_height)
    })
    return(c(heatmaps, violins))
}

mk_cluster_counts_data  <- function(liger_meta) {
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

# format cluster dge values.
# anno_points will rearrange columns but ONLY if the matrix ordering matches the unclustered sorting order of the main heatmap matrix.
# therefore, reorder by dge_hclust$labels, which shuld have the lexicographical order
mk_cluster_limma_data <- function(dge_data, dge_hclust, limma_tb) {
    subclusters <- dge_data %>% summarize(ct_subcluster = as.character(unique(ct_subcluster)))
    dge_order <- dge_hclust$labels
    limma_filter <- limma_tb %>%
        select(-data, -matrix, -limma, -limma_bayes) %>%
        unnest(limma_pval) %>%
        right_join(subclusters, by = c("ct_subcluster"))

    # format single dx v. control p values.
    limma_single_dx <- pivot_longer(limma_filter, cols = matches("^pvaldx.*Control$"), names_to = "type", values_to = "value") %>%
        select(ct_subcluster, type, value) %>%
        mutate(value = -log10(value))

    # format single dx v. all p values.
    limma_all_ct_dx <- pivot_longer(limma_filter, cols = matches("^pvaldx.*/ 3$"), names_to = "type", values_to = "value") %>%
        select(ct_subcluster, type, value) %>%
        mutate(value = -log10(value))

    # format as matrix form; reorder rows according to dge hclust
    limma_single_dx_mat <- pivot_matrix(limma_single_dx, "type", "value", "ct_subcluster")
    limma_all_ct_mat <- pivot_matrix(limma_all_ct_dx, "type", "value", "ct_subcluster")
    limma_single_dx_mat <- limma_single_dx_mat[dge_order, ]
    limma_all_ct_mat <- limma_all_ct_mat[dge_order, ]

    return(list(limma_single_dx_mat = limma_single_dx_mat, limma_all_ct_mat = limma_all_ct_mat))
}

pal_stallion <- c("1"="#D51F26","2"="#272E6A","3"="#208A42","4"="#89288F","5"="#F47D2B", "6"="#FEE500","7"="#8A9FD1","8"="#C06CAB","19"="#E6C2DC",
                "10"="#90D5E4", "11"="#89C75F","12"="#F37B7D","13"="#9983BD","14"="#D24B27","15"="#3BBCA8", "16"="#6E4B9E","17"="#0C727C", "18"="#7E1416","9"="#D8A767","20"="#3D3D3D")

# make annotations based on cell counts.
# adding a legend using annotation_legend_param is not supported for anno_barplot -- return as list to plot later down the line
mk_cluster_counts_annot <- function(count_data) {
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

# make dotplot showing limma results.
mk_cluster_limma_annot <- function(limma_data) {
    limma_single_dx_mat <- limma_data$limma_single_dx_mat
    limma_all_ct_mat <- limma_data$limma_all_ct_mat

    # pick from pal_stallion palette.
    single_dx_colors <- pal_stallion[1:ncol(limma_single_dx_mat)]
    all_ct_colors <- pal_stallion[1:ncol(limma_all_ct_mat)]
    
    ylim_common <- c(0, max(limma_single_dx_mat, limma_all_ct_mat, na.rm = TRUE))

    return(
        list(
            annotations = list(
                single_dx = columnAnnotation(
                    single_dx = anno_points(limma_single_dx_mat, which = "column", gp = gpar(col = single_dx_colors), ylim = ylim_common),
                    annotation_label = "limma\ndx v. ctl",
                    annotation_name_rot = 90,
                    annotation_height = unit(20, "mm")
                ),
                all_ct = columnAnnotation(
                    all_ct = anno_points(limma_all_ct_mat, which = "column", gp = gpar(col = all_ct_colors), ylim = ylim_common),
                    annotation_label = "limma\ndx v. all",
                    annotation_name_rot = 90,
                    annotation_height = unit(20, "mm")
                )
            ),
            legends = list(
                single_dx = Legend(labels = colnames(limma_single_dx_mat), title = "-log10 limma pval single dx v. control", legend_gp = gpar(fill = single_dx_colors)),
                all_ct = Legend(labels = colnames(limma_all_ct_mat), title = "-log10 limma pval single dx v. all", legend_gp = gpar(fill = all_ct_colors))
            )
        )
    )

}

mk_nd_heatmap <- function(nd_data) {
    nd_mtx <- nd_data %>% reduce(., ~inner_join(.x, .y, by = "ct_subcluster")) %>%
        column_to_rownames("ct_subcluster") %>%
        as.matrix() %>%
        t()

    colormap <- colorRamp2(
        breaks = c(min(nd_mtx, na.rm = TRUE), 0, max(nd_mtx, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )

    Heatmap(
        nd_mtx,
        row_title = " (b) ",
        col = colormap,
        height = unit(20, "mm"),
        cluster_columns = FALSE,
        cluster_rows = FALSE
    )
}

mk_combo_heatmap <- function(...) {
    cr <- list(...)
    hmap <- cr$dge_base_heatmap
    hmap <- hmap %v% cr$cluster_nd_heatmap
    for (annot in cr$gene_expr_annot) {
        hmap <- hmap %v% annot
    }
    for (annot in cr$cluster_counts_annot$annotations) {
        hmap <- hmap %v% annot
    }
    for (annot in cr$cluster_limma_annot$annotations) {
        hmap <- hmap %v% annot
    }
    # names are dropped after adding annotations. recover by manually adding text annotation.
    # col_labels are reordered during heatmap draw -- don't reorder manually
    base_col_labels <- cr$dge_hclust$labels
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

if (!interactive()) {
    main()
}
