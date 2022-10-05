# Complexheatmap showing clustering of liger clusters with:
# 1. violin plot of expression values of selected genes.
# 2. counts of genes / cells per cluster.

set.seed(0)
options(bitmapType='png')
liblist <- c("Seurat", "tidyverse", "readxl", "ComplexHeatmap", "circlize", "RColorBrewer", "scales", "WGCNA", "cowplot", "broom.mixed", "patchwork")
l <- lapply(liblist, function(x) {suppressPackageStartupMessages(require(x, character.only = TRUE, quietly = TRUE))})
in_cluster_dge_wk <- "../../analysis/seurat_lchen/liger_subcluster_enrichment_dge/subcluster_wk.rds"
in_cluster_dx_dge_wk <- "../../analysis/seurat_lchen/liger_subcluster_lme/subcluster_wk.rds"
in_limma <- "../../analysis/seurat_lchen/subcluster_composition/limma/subcluster_composition_dx_pmi.rds"
in_seurat_liger <- "../../analysis/seurat_lchen/liger_subcluster_metadata.rds"
in_seurat_rds <- "../../analysis/pci_import/pci_seurat.rds"

clusters_exclude_file <- "../../resources/subclusters_removed_byQC_final.xlsx"
#excitatory_markers <- "../../resources/excitatory_layers_20210318.csv"
marker_file <- "../../resources/CelltypeMarkers_biologicalGroupings.xlsx"
in_bulk_meta <- "../../resources/individual_nd.rds"

out_path_base <- "../../analysis/seurat_lchen/liger_subcluster_hier/heatmap"
batchtools <- file.path(out_path_base, "batchtools")
options(future.globals.maxSize = Inf)

main <- function() {
    dir.create(out_path_base, recursive = TRUE, showWarnings = FALSE)

    writeLines("read")
    seurat_obj <- readRDS(in_seurat_rds)
    cluster_wk_in <- readRDS(in_cluster_dge_wk)
    cluster_wk_dx_in <- readRDS(in_cluster_dx_dge_wk)
    liger_meta_in <- readRDS(in_seurat_liger)
    bulk_meta <- readRDS(in_bulk_meta)
    excludes <- read_xlsx(clusters_exclude_file)
    limma_tb <- readRDS(in_limma)
    marker_tb <- read_xlsx(marker_file, sheet = "combined_markers") %>%
        mutate(group = as.factor(replace_na(group, "markers"))) %>%
        group_by(cluster_cell_type) %>%
        group_nest(.key = "markers")

    nd_tb <- bulk_meta %>%
        mutate(Autopsy.ID = paste0("P", Autopsy.ID)) %>%
        group_by(Autopsy.ID, type) %>%
        slice_head(n = 1) %>%
        select(Autopsy.ID, type, score) %>%
        mutate(Autopsy.ID = factor(Autopsy.ID), type = factor(type))

    cluster_dge_wk <- cluster_wk_in %>%
        filter(!ct_subcluster %in% excludes$ct_subcluster) %>%
        filter(!is.na(broom_join))
    
    cluster_ct_group <- cluster_dge_wk %>%
        group_by(cluster_cell_type) %>%
        group_nest(.keep = T) %>%
        filter(!cluster_cell_type %in% c("t_cell"))

    cluster_ct_group <- mutate(cluster_ct_group, liger_meta = map(data, function(data) {
        filter(liger_meta_in, ct_subcluster %in% data$ct_subcluster)
    }))

    writeLines("get enrichment broom data")
    cluster_ct_group <- mutate(cluster_ct_group,
        dge_data = map(data, ~tryCatch(fmt_cluster_dge(.), error = function(x) {print(x); return(NA)}))
    ) %>%
    filter(!is.na(dge_data))
    
    cluster_ct_group <- cluster_ct_group %>%
        select(-any_of("markers")) %>%
        left_join(marker_tb, by = "cluster_cell_type")


    writeLines("filter top 100 variable genes per subcluster")
    cluster_ct_group <- mutate(cluster_ct_group, dge_plot_tb = map(dge_data, filter_var_genes, filter_genes_all_regions = TRUE))
    writeLines("hclust per seurat celltype on variable genes")
    cluster_ct_group <- mutate(cluster_ct_group, dge_hclust = map(dge_plot_tb, mk_dge_hclust))
    writeLines("get gene expression per marker")
    cluster_ct_group <- mutate(cluster_ct_group, gene_expr_data = pmap(list(liger_meta, markers), mk_gene_data, seurat_obj))
    writeLines("format enrichment dge per cluster")
    cluster_ct_group <- mutate(cluster_ct_group, dge_dx_data = pmap(list(data, markers), fmt_cluster_dx_dge, cluster_wk_dx_in))
    writeLines("get cell counts per subcluster")
    cluster_ct_group <- mutate(cluster_ct_group, cluster_counts_data = map(liger_meta, mk_cluster_counts_data))
    writeLines("format betas as matrix")
    cluster_ct_group <- mutate(cluster_ct_group, dge_marker_data = pmap(list(dge_data, dge_hclust, markers), mk_beta_marker_data))
    writeLines("format limma subcluster data")
    cluster_ct_group <- mutate(cluster_ct_group, cluster_limma_data = pmap(list(dge_data, dge_hclust), mk_cluster_limma_data, limma_tb))
    writeLines("run bicor of nd data against cluster cell counts")
    cluster_ct_group <- mutate(cluster_ct_group, nd_correlation = map(liger_meta, mk_nd_data, nd_tb))
    
    writeLines("plot format")
    cluster_ct_group <- mutate(cluster_ct_group, dge_base_heatmap = pmap(list(dge_plot_tb, dge_hclust), mk_beta_heatmap))
    cluster_ct_group <- mutate(cluster_ct_group, gene_expr_annot = map(gene_expr_data, mk_gene_annot))
    cluster_ct_group <- mutate(cluster_ct_group, dge_dx_annot = pmap(list(dge_dx_data, dge_data, dge_hclust), mk_dge_dx_annot))
    cluster_ct_group <- mutate(cluster_ct_group, dge_marker_heatmap = pmap(list(dge_marker_data, dge_hclust, markers), mk_beta_marker_heatmaps))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_counts_annot = map(cluster_counts_data, mk_cluster_counts_annot))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_limma_annot = map(cluster_limma_data, mk_cluster_limma_annot))
    cluster_ct_group <- mutate(cluster_ct_group, cluster_nd_heatmap = map(nd_correlation, mk_nd_heatmap))
    cluster_ct_group <- mutate(cluster_ct_group, dge_combo_heatmap = pmap(cluster_ct_group, mk_combo_heatmap))

    writeLines("write pdf")
    pwalk(cluster_ct_group, function(...) {
        cr <- list(...)
        heatmap_path <- file.path(out_path_base, str_glue("dge_enrichment_{cr$cluster_cell_type}_var_beta_heatmap.pdf"))
        plot_ts <- str_glue("{system('md5sum liger_subcluster_complexheatmap.R', intern = TRUE)} {date()}")

        tmp_pdf <- function(w, h) { pdf(tempfile(), w, h) }

        if (nrow(cr$dge_marker_data) > 1) {
            height <- unit((0.5 * nrow(cr$dge_marker_data) + nrow(cr$markers)) + 2, "cm")
        } else {
            height <- unit(10, "cm")
        }
        heatmap_gtree <- grid.grabExpr(
            draw(cr$dge_combo_heatmap, heatmap_legend_list = c(cr$cluster_counts_annot$legends, cr$cluster_limma_annot$legends)),
            wrap = TRUE,
            device = tmp_pdf,
            width = unit(10, "cm"),
            height = height
        )

        print(heatmap_path)
        pdf(heatmap_path, width = unit(15, "cm"), height = height)
        print(wrap_plots(heatmap_gtree) + plot_annotation(title = cr$cluster_cell_type, subtitle = plot_ts))
        graphics.off()
    })

}

# data processing funcs

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

fmt_cluster_dx_dge <- function(cluster_data, markers, cluster_dx_wk) {
    dx <- c("AD", "bvFTD", "PSP-S")
    clinical_dx_beta_cols <- str_glue("clinical_dx{dx}.estimate")
    clinical_dx_pval_cols <- str_glue("clinical_dx{dx}.p.value")
    clinical_dx_fdr_cols <- str_glue("clinical_dx{dx}.p.value.adj")
    clinical_dx_z_cols <- str_glue("clinical_dx{dx}.statistic")
    filter_args <- tibble(dx, clinical_dx_beta_cols, clinical_dx_pval_cols, clinical_dx_fdr_cols, clinical_dx_z_cols)
   
    filtered_gene_tests <- cluster_data %>%
        select(ct_subcluster) %>%
        left_join(cluster_dx_wk, by = c("ct_subcluster")) %>%
        mutate(ct_subcluster = fct_drop(ct_subcluster)) %>%
        full_join(., filter_args, by = character())

    filter_call <- function(...) {
        cr <- list(...)
        beta_col <- cr$clinical_dx_beta_cols
        pval_col <- cr$clinical_dx_pval_cols
        fdr_col <- cr$clinical_dx_fdr_cols
        z_col <- cr$clinical_dx_z_cols
        writeLines(str_glue("{cr$ct_cluster}:  {beta_col} {pval_col} {fdr_col} {z_col}"))
        if (!is.na(cr$broom_join) && nrow(cr$broom_join) > 1 && all(c(beta_col, fdr_col, z_col) %in% colnames(cr$broom_join))) {
            broom_filter <- cr$broom_join %>%
                dplyr::filter(gene %in% markers$gene) %>%
                rename(beta = .data[[beta_col]], pval = .data[[pval_col]], fdr = .data[[fdr_col]], z = .data[[z_col]]) %>%
                select(gene, beta, pval, fdr, z)
        }
    }
    filtered_gene_tb <- filtered_gene_tests %>%
        mutate(lme_marker_estimates = pmap(., filter_call)) %>%
        select(region, dx, ct_subcluster, lme_marker_estimates) %>%
        unnest(lme_marker_estimates) %>%
        complete(ct_subcluster, dx)
    return(filtered_gene_tb) 
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
    dge_beta_matrix <- dge_beta_matrix[, beta_hclust$labels]
    colormap <- colorRamp2(
        breaks = c(min(dge_beta_matrix, na.rm = TRUE), 0, max(dge_beta_matrix, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )
    heat <- Heatmap(
        dge_beta_matrix,
        name = "gene beta hclust",
        col = colormap,
        na_col = "grey75",
        cluster_columns = beta_hclust,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        height = unit(20, "mm")
    )
    return(heat)
}

mk_beta_marker_data <- function(dge_data, beta_hclust, annot_gene_list) {
    # dge data may not have results for all subclusters; using complete() to match beta_hclust subclusters
    gene_dge_list <- dge_data %>%
        filter(gene %in% annot_gene_list$gene) %>%
        group_by(gene) %>%
        complete(ct_subcluster = beta_hclust$labels)
    dge_beta_matrix <- pivot_matrix(gene_dge_list, "ct_subcluster", "beta", "gene")
    dropped_genes <- setdiff(annot_gene_list$gene, rownames(dge_beta_matrix))
    writeLines(str_glue("{dge_data$data[[1]]$cluster_cell_type[[1]]} dropped genes: ", paste0(dropped_genes, collapse = " ")))
    if (nrow(dge_beta_matrix) == 0) return(matrix())
    stopifnot(all(colnames(dge_beta_matrix) %in% beta_hclust$labels))
    dge_beta_matrix <- dge_beta_matrix[, beta_hclust$labels]
    return(dge_beta_matrix)
}

# - Subset seurat object to cells in dge_plot_tb / genes in annot_gene_list.
# - Regress out variables for display purposes using negative binomial model.
# - For each gene, get expression matrix and convert to tibble with cols expr / umi.
#   - join with cell types' metadata and split along ct_subcluster to get a list of vectors containing each subclusters' expression values.
# Returns a list of lists:
# - names(x) is gene names in annot_gene_list.
# - names(x[[*]]) is subclusters in liger_meta.
# - x[[*]] is a vector of gene expression values for each subcluster.
mk_gene_data <- function(liger_meta, annot_gene_list, sobj) {
    feature_list <- intersect(annot_gene_list$gene, rownames(GetAssayData(sobj, slot = "scale.data")))
    subset_sobj <- subset(sobj, cells = liger_meta$UMI, features = feature_list)
    #sobj_scaled <- ScaleData(subset_sobj, model.use = "negbinom", vars.to.regress = c("pmi", "age", "sex", "number_umi", "percent_mito"), assay = "RNA")
    sobj_mat <- GetAssayData(subset_sobj, slot = "scale.data")

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

symnum_signif <- function(pval, pval.fdr) {
    ifelse(
        pval.fdr < 0.1,
        symnum(pval.fdr, cutpoints = c(0, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", " ")),
        ifelse(
            pval < 0.05,
            "#",
            " "
        )
    )
}

# Run bicor of nd values against cluster cell counts.
mk_nd_data <- function(liger_meta, nd_tb) {
    library_celltype_counts_full <- liger_meta %>%
        filter(cluster_cell_type == cell_type) %>%
        mutate(ct_subcluster = fct_drop(ct_subcluster)) %>%
        group_by(library_id, ct_subcluster) %>%
        summarize(
            autopsy_id = unique(autopsy_id),
            clinical_dx = unique(clinical_dx),
            age = unique(age),
            pmi = unique(pmi),
            sex = unique(sex),
            log_pmi = log(unique(pmi)),
            mean_percent_mito = mean(percent_mito),
            median_genes = median(number_genes),
            cluster_ct = n(),
            .groups = "drop"
        )

    # bicor / cor.test on cell counts per subcluster.
    meta_nd_tb <- left_join(library_celltype_counts_full, nd_tb, by = c("autopsy_id" = "Autopsy.ID")) %>%
        filter(!is.na(type)) %>%
        group_by(ct_subcluster, type) %>%
        summarize(
            bicor = bicor(cluster_ct, score, use = "pairwise.complete.obs", pearsonFallback = "none")[, 1],
            cor.test = tidy(cor.test(cluster_ct, score)),
            cor.pval = cor.test$p.value,
            cor.fdr = p.adjust(cor.pval),
            cor.signif = symnum_signif(cor.pval, cor.fdr),
            .groups = "drop"
        )
    # fill in subclusters that didn't have a sample in nd_tb
    meta_nd_tb <- meta_nd_tb %>% complete(expand(meta_nd_tb, ct_subcluster, type))


    return(meta_nd_tb)
}

# heatmap formatting funcs

mk_beta_marker_heatmaps <- function(marker_data, beta_hclust, marker_tb, row_hclust = FALSE) {
    if (!is.null(marker_tb)) {
        marker_tb %>%
            group_by(group) %>%
            group_map(function(.x, .y, ...) {
                genes <- intersect(rownames(marker_data), .x$gene)
                gene_order <- .x$gene[.x$gene %in% genes]
                marker_ord_data <- marker_data[gene_order,,drop = F]
                if (is.matrix(marker_ord_data) && nrow(marker_ord_data) > 0) {
                    colormap <- colorRamp2(
                        breaks = c(-2, 0, 2),
                        colors = c(muted("blue"), "white", muted("red"))
                    )
                    height <- unit(nrow(marker_ord_data) * 10, "mm")
                    heat <- Heatmap(
                        marker_ord_data,
                        name = "gene beta",
                        col = colormap,
                        na_col = "grey75",
                        cluster_columns = beta_hclust,
                        cluster_rows = FALSE,
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        height = height
                    )
                    return(heat)
                } else {
                    return(NULL)
                }
            })
    }

}

mk_gene_annot <- function(gene_data) {
    annotation_height <- unit(10, "mm")
    heatmaps <- imap(gene_data, function(subcluster_expr_vec_list, gene_name){
        columnAnnotation(gene = anno_density(subcluster_expr_vec_list, type = "heatmap"), annotation_label = gene_name, annotation_height = annotation_height)
    })
    violins <- imap(gene_data, function(subcluster_expr_vec_list, gene_name){
        columnAnnotation(gene = anno_density(subcluster_expr_vec_list, type = "violin"), annotation_label = gene_name, annotation_height = annotation_height)
    })
    #return(c(heatmaps, violins))
    return(heatmaps)
}

pal_cts_stallion <- c("AD"="#D51F26","bvFTD"="#272E6A","Control"="#89288F","PSP-S"="#208A42") # override counts color for control

mk_dge_dx_annot <- function(dge_dx_data, dge_data, dge_hclust) {
    subclusters <- dge_data %>% summarize(ct_subcluster = as.character(unique(ct_subcluster)))
    dge_order <- dge_hclust$labels
    dx_order <- unique(dge_dx_data$dx)
    pal <- pal_cts_stallion[dx_order]
    annotation_height <- unit(10, "mm")

    dge_dx_mats <- dge_dx_data %>%
        mutate(ct_subcluster = fct_expand(ct_subcluster, dge_order), dx = fct_expand(dx, dx_order)) %>%
        group_split(gene) %>%
        map(function(gene_tb) {
            gene_mat <- gene_tb %>%
                complete(ct_subcluster, dx) %>%
                pivot_matrix("ct_subcluster", "beta", "dx")
        }) %>%
        setNames(unique(dge_dx_data$gene))


    annots <- imap(dge_dx_mats, function(dge_mat, name) {
        if (!is.na(dge_mat)) {
            columnAnnotation(dge_dx = anno_points(t(dge_mat), gp = gpar(col = pal)), annotation_label = name, annotation_height = annotation_height)
        }
    })
    return(annots)
}

# compute cluster 
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

    # format single dx v. control beta.
    limma_single_dx_beta <- limma_tb %>%
        select(-data, -matrix, -limma, -limma_bayes) %>%
        unnest(limma_beta) %>%
        right_join(subclusters, by = c("ct_subcluster")) %>%
        pivot_longer(cols = matches("^betadx.*Control$"), names_to = "type", values_to = "value") %>%
        select(ct_subcluster, type, value)

    # format single dx v. all p values.
    limma_all_dx_beta <- limma_tb %>%
        select(-data, -matrix, -limma, -limma_bayes) %>%
        unnest(limma_beta) %>%
        right_join(subclusters, by = c("ct_subcluster")) %>%
        pivot_longer(cols = matches("^betadx.*/ 3$"), names_to = "type", values_to = "value") %>%
        select(ct_subcluster, type, value)

    # format as matrix form; reorder rows according to dge hclust
    # multiply by sign() of estimate 

    limma_single_dx_est_mat <- pivot_matrix(limma_single_dx_beta, "type", "value", "ct_subcluster")
    limma_all_ct_est_mat <- pivot_matrix(limma_all_dx_beta, "type", "value", "ct_subcluster")

    limma_single_dx_mat <- pivot_matrix(limma_single_dx, "type", "value", "ct_subcluster")
    limma_all_ct_mat <- pivot_matrix(limma_all_ct_dx, "type", "value", "ct_subcluster")

    limma_single_dx_mat <- limma_single_dx_mat * sign(limma_single_dx_est_mat)
    limma_all_ct_mat <- limma_all_ct_mat * sign(limma_all_ct_est_mat)

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
    dx_annot_colors <- pal_cts_stallion[1:length(dx_names)]
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
    
    ylim_common <- c(min(limma_single_dx_mat, limma_all_ct_mat, na.rm = TRUE), max(limma_single_dx_mat, limma_all_ct_mat, na.rm = TRUE))

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
    nd_mtx <- pivot_matrix(nd_data, "ct_subcluster", "bicor", "type")
    heatmap_text_matrix <- pivot_matrix(nd_data, "ct_subcluster", "cor.signif", "type")

    colormap <- colorRamp2(
        breaks = c(min(nd_mtx, na.rm = TRUE), 0, max(nd_mtx, na.rm = TRUE)),
        colors = c(muted("blue"), "white", muted("red"))
    )
    
    plot_text_label <- function(j, i, x, y, w, h, col) {
        label <- heatmap_text_matrix[i, j]
        if (!is.na(label)) {
            grid.text(label, x, y)
        }
    }

    Heatmap(
        nd_mtx,
        name = "library count bicor",
        row_title = " (b) ",
        cell_fun = plot_text_label,
        col = colormap,
        height = unit(20, "mm"),
        cluster_columns = FALSE,
        cluster_rows = FALSE
    )
}

mk_combo_heatmap <- function(...) {
    cr <- list(...)
    print(cr$cluster_cell_type)
    hmap <- cr$dge_base_heatmap
    hmap <- hmap %v% cr$cluster_nd_heatmap
    for (annot in cr$dge_dx_annot) {
        hmap <- hmap %v% annot
    }
    for (annot_heatmap in cr$dge_marker_heatmap) {
        hmap <- hmap %v% annot_heatmap
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

