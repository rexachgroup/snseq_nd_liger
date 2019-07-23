# Damon Polioudakis
# 2018-10-18
# Plots of seurat analysis

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)

require(methods)
require(Seurat)
require(Matrix)
require(cowplot)
require(loomR)
require(ggcorrplot)
require(viridis)
require(RColorBrewer)
require(tidyverse)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
sessionInfo()
# require(xlsx)

## inputs
in_seurat <- paste0(
  "/u/flashscratch/d/dpolioud/seurat_subcluster/20190528/cca10_lib/"
  , "seurat_subcluster_astro.rdat")
in_seurat_raw <- paste0(
  "/u/flashscratch/d/dpolioud/seurat/20190427/"
  , "p1_p2_p3_p4_raw.rdat")
  # , "/p1_p2_raw_test.rdat")
marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
astro_markers_tb <- read_csv(
  "../resources/AstrocyteMarkers_Brie_20190110.csv")

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- "seurat_plots.R"
# graph_subtitle <- "P1 P2 P3 P4 astrocytes batch correct with CCA (10 CCs) on library"
graph_subtitle <- "P1 P2 P3 P4"

## outputs
out_graph <- paste0("../analysis/seurat/", date
  , "/graphs/p1_p2_p3_p4/seurat_")
out_table <- paste0("../analysis/seurat/", date
  , "/tables/p1_p2_p3_p4/seurat_")
# out_graph <- paste0("../analysis/seurat_subcluster/", date
#   , "/graphs/p1_p2_p3_p4/cca10_lib/seurat_")
# out_table <- paste0("../analysis/seurat_subcluster/", date
#   , "/tables/p1_p2_p3_p4/cca10_lib/seurat_")

# make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
################################################################################

### functions

main_function <- function(){
  print("main_function")
  load(in_seurat)
  # for when Seurat v2.3 was used for running CCA (not yet in v3.0)
  # p1_p2_so <- UpdateSeuratObject(subset_so)
  p1_p2_so <- p1_p2_p3_p4_so
  rm(p1_p2_p3_p4_so)
  # cluster_enriched_de_tb <- cluster_enriched_tb
  p1_p2_so$cluster_pc1to100_res06 <- p1_p2_so[["RNA_snn_res.0.6"]]
  p1_p2_so$cell_ids <- p1_p2_so[[]] %>% rownames
  # plot_qc_metrics(seurat_obj = p1_p2_so, in_seurat_raw = in_seurat_raw)
  plot_pca(seurat_obj = p1_p2_so)
  plot_correlation_matrix_of_pcs_and_metadata(seurat_obj = p1_p2_so)
  # plot_tsne_clustering_pcs_test()
  # plot_tsne_clustering_resolution_test(seurat_obj = p1_p2_so)
  plot_dim_reduction_colored_by_metadata(
    seurat_obj = p1_p2_so, reduction = "tsne")
  plot_dim_reduction_colored_by_metadata(
    seurat_obj = p1_p2_so, reduction = "umap")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = p1_p2_so, reduction = "tsne")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = p1_p2_so, reduction = "umap")
  plot_tsne_colored_by_pc_score(seurat_obj = p1_p2_so)
  plot_cluster_metrics(seurat_obj = p1_p2_so)
  plot_marker_expression_dim_reduction(
    seurat_obj = p1_p2_so, reduction = "tsne")
  plot_marker_expression_dim_reduction(
    seurat_obj = p1_p2_so, reduction = "umap")
  # make_table_of_percent_of_cell_types()
  plot_marker_astrocyte_expression_dim_reduction(
    seurat_obj = subset_so, reduction = "tsne")
  plot_marker_astrocyte_expression_dim_reduction(
    seurat_obj = subset_so, reduction = "umap")
  # plot_genes_of_interest_expression_tsnes(seurat_obj = p1_p2_so)
  # plot_cluster_enriched_genes_heatmap(
    # seurat_obj = p1_p2_so, cluster_enriched_de_tb = cluster_enriched_de_tb)
  plot_cluster_top_expressed_genes_heatmap(
    seurat_obj = p1_p2_so, top_expressed_genes_tb = top_expressed_genes_tb)
  # plot_cluster_enriched_genes_violin_plots(
  #   seurat_obj = p1_p2_so, cluster_enriched_de_tb = cluster_enriched_de_tb)
}

main_function_subclustered_astrocytes <- function(){
  print("main_function_subclustered_astrocytes")
  load(in_seurat)
  # for when Seurat v2.3 was used for running CCA (not yet in v3.0)
  subset_so <- UpdateSeuratObject(subset_so)
  # cluster_enriched_de_tb <- cluster_enriched_tb
  subset_so$subclustering <- subset_so[["RNA_snn_res.0.4"]]
  subset_so$cell_ids <- subset_so[[]] %>% rownames
  # plot_qc_metrics(seurat_obj = subset_so, in_seurat_raw = in_seurat_raw)
  # plot_pca(seurat_obj = subset_so)
  # plot_correlation_matrix_of_pcs_and_metadata(seurat_obj = subset_so)
  # plot_tsne_clustering_pcs_test()
  # plot_tsne_clustering_resolution_test(seurat_obj = subset_so)
  plot_dim_reduction_colored_by_metadata(
    seurat_obj = subset_so, reduction = "tsne")
  plot_dim_reduction_colored_by_metadata(
    seurat_obj = subset_so, reduction = "umap")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = subset_so, cluster_col_name = "cluster_pc1to100_res06"
    , reduction = "tsne")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = subset_so
    , cluster_col_name = "subclustering", reduction = "tsne")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = subset_so, cluster_col_name = "cluster_pc1to100_res06"
    , reduction = "umap")
  plot_dim_reduction_colored_by_cluster(
    seurat_obj = subset_so
    , cluster_col_name = "subclustering", reduction = "umap")
  # plot_tsne_colored_by_pc_score(seurat_obj = subset_so)
  plot_cluster_metrics(
    seurat_obj = subset_so, cluster_col_name = "subclustering")
  plot_marker_expression_dim_reduction(
    seurat_obj = subset_so, cluster_col_name = "subclustering"
    , reduction = "tsne")
  plot_marker_expression_dim_reduction(
    seurat_obj = subset_so, cluster_col_name = "subclustering"
    , reduction = "umap")
  # make_table_of_percent_of_cell_types()
  plot_marker_astrocyte_expression_dim_reduction(
    seurat_obj = subset_so, cluster_col_name = "subclustering"
    , reduction = "tsne")
  plot_marker_astrocyte_expression_dim_reduction(
    seurat_obj = subset_so, cluster_col_name = "subclustering"
    , reduction = "umap")
  # plot_genes_of_interest_expression_tsnes(seurat_obj = subset_so)
  # plot_cluster_enriched_genes_heatmap(
    # seurat_obj = subset_so, cluster_enriched_de_tb = cluster_enriched_de_tb)
  plot_cluster_top_expressed_genes_heatmap(
    seurat_obj = subset_so, top_expressed_genes_tb = top_expressed_genes_tb)
  # plot_cluster_enriched_genes_violin_plots(
  #   seurat_obj = subset_so, cluster_enriched_de_tb = cluster_enriched_de_tb)
}
################################################################################

### plot QC metrics (histograms, density plots)

plot_qc_metrics <- function(seurat_obj, in_seurat_raw){

  print("plot_qc_metrics")

  load(in_seurat_raw)
  p1_p2_p3_p4_raw_so$cell_ids <- p1_p2_p3_p4_raw_so[[]] %>% rownames

  ## Plots

  metadata <- c("cell_ids", "number_genes", "number_umi", "percent_mito"
    , "library_id")
    # , "mean_gc_content", "mean_cds_length")
  dat_tb <- p1_p2_p3_p4_raw_so[[metadata]] %>%
    # rename(. ,replace = c("sample_number" = "sample_number")) %>%
    mutate(data_type = "raw data") %>%
      bind_rows(.
        , seurat_obj[[metadata]] %>%
          mutate(data_type = "filtered data")
        ) %>%
      mutate(data_type = factor(
        data_type, levels = c("raw data", "filtered data"))) %>%
      as_tibble() %>%
      mutate(mean_gc_content = 1) %>%
      mutate(mean_cds_length = 1)
      # gather(key = metric, value = value, -data_type, -cell_ids)

  # violin boxplots
  means_tb <- dat_tb %>%
  gather(variable, value, -data_type, -cell_ids, -library_id) %>%
    group_by(data_type, variable) %>%
    summarize(
      value_mean = round(mean(value, na.rm = TRUE), 1)
      , value_median = round(median(value, na.rm = TRUE),1)
    ) %>%
    mutate(label = paste0(
      "mean:\n", value_mean, "\nmedian:\n", value_median))
  dat_tb %>%
    gather(variable, value, -data_type, -cell_ids, -library_id) %>%
    group_by(data_type, variable) %>%
    ggplot(aes(x = data_type, y = value, fill = data_type)) +
      facet_wrap(~variable, scales = "free", ncol = 5) +
      geom_violin(fill = "lightgrey", color = "lightgrey") +
      geom_boxplot(width = 0.2, outlier.size = 0.1) +
      theme(
        legend.position = "none"
        , axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data = means_tb, aes(
        x = data_type, y = Inf, label = label), vjust = 1) +
      ggtitle(make_plot_title("QC metrics"))
    ggsave(paste0(out_graph, "qc_metrics_violin_boxplots.pdf")
      , width = 9, height = 5)

  plot_grid_wrapper(
    list(
      # number genes
      # histograms
      ggplot(dat_tb, aes(x = number_genes)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_genes, color = data_type)) +
        geom_density()
      # number umi
      # histograms
      , ggplot(dat_tb, aes(x = number_umi)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = number_umi, color = data_type)) +
        geom_density() +
        ggtitle("Number UMI (>0 counts). Binwidth: 50")
      # percent MT
      # histograms
      , ggplot(dat_tb, aes(x = percent_mito)) +
        geom_histogram(binwidth = 1) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = percent_mito, color = data_type)) +
        geom_density() +
        ggtitle("Percentage MT. Binwidth: 1")

      # mean CDS length
      # histograms
      , ggplot(dat_tb, aes(x = mean_cds_length)) +
        geom_histogram(binwidth = 50) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = mean_cds_length, color = data_type)) +
        geom_density() +
        ggtitle("Mean CDS length. Binwidth: 50")
      # mean gc content
      # histograms
      , ggplot(dat_tb, aes(x = mean_gc_content)) +
        geom_histogram(binwidth = 1) +
        facet_wrap(~data_type, ncol = 2)
      # density plots
      , ggplot(dat_tb, aes(x = mean_gc_content, color = data_type)) +
        geom_density() +
        ggtitle("Percentage GC. Binwidth: 1")
      )
    , rel_height = 0.1
    , title = (make_plot_title("Histograms and density plots of QC metrics"
      ))
  )
  ggsave(paste0(out_graph, "qc_metrics_histogram_density.png")
    , width = 13, height = 16)

  # number_genes
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = number_genes, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        coord_cartesian(xlim = c(0,6000)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_genes, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,6000)) +
          facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(paste0(
      "Histograms and density plots: number of genes detected per cell by sample"
      , "\nx axis limits 0-6,000"
    )))
  )
  ggsave(paste0(out_graph, "qc_number_genes_by_sample_histogram_density.png")
    , width = 22, height = 6)

  # number_umi
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = number_umi, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        coord_cartesian(xlim = c(0,10000)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Number genes (>0 counts). Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = number_umi, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,10000)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(paste0(
      "Histograms and density plots: number of UMI per cell by sample"
      , "\nx axis limits 0-10,000"
    )))
  )
  ggsave(paste0(out_graph, "qc_number_umi_by_sample_histogram_density.png")
    , width = 22, height = 5)

  # percent_mito
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = percent_mito, fill = library_id)) +
        geom_histogram(binwidth = 1) +
        facet_grid(data_type~library_id) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Percentage MT. Binwidth: 1")
      # density plots
      , ggplot(dat_tb, aes(x = percent_mito, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        coord_cartesian(xlim = c(0,25)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: percent MT per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_percent_mt_by_sample_histogram_density.png")
    , width = 22, height = 6)

  # mean cds length
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = mean_cds_length, fill = library_id)) +
        geom_histogram(binwidth = 50) +
        facet_grid(data_type~library_id) +
        # scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Mean CDS length. Binwidth: 50")
      # density plots
      , ggplot(dat_tb, aes(x = mean_cds_length, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: mean CDS length per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_mean_cds_length_by_sample_histogram_density.png")
    , width = 22, height = 5)

  # mean gc content
  plot_grid_wrapper(
    list(
      # histograms
      ggplot(dat_tb, aes(x = mean_gc_content, fill = library_id)) +
        geom_histogram(binwidth = 1) +
        facet_grid(data_type~library_id) +
        theme(
          legend.position = "none"
          , axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Mean GC content. Binwidth: 1")
      # density plots
      , ggplot(dat_tb, aes(x = mean_gc_content, color = library_id)) +
        geom_density(size = 0.5, aes(linetype = library_id)) +
        scale_linetype_manual(values = rep(1:12,100)) +
        facet_wrap(~data_type, ncol = 1)
      )
    , rel_height = 0.2, rel_widths = c(1, 0.5), align = 'h', axis = 't'
    , title = (make_plot_title(
      "Histograms and density plots: mean gc content per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_mean_gc_by_sample_histogram_density.png")
    , width = 22, height = 6)

}
################################################################################

### PCA plots

plot_pca <- function(seurat_obj){

  print("plot_pca")

  # Plot genes with highest PC loadings
  pc_df <- Loadings(object = seurat_obj, reduction = "pca")[,1:8] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble %>%
    gather(., key = "pc", value = "pc_loading", -gene) %>%
    mutate(pc = gsub("V", "", pc)) %>%
    group_by(pc) %>%
    filter(! is.na(pc_loading)) %>%
    arrange(pc_loading) %>%
    filter(row_number() %in% c(1:15,(n() - 14):n())) %>%
    ungroup() %>%
    # filter(pc == "PC_1")
    as.data.frame()
  gg_l <- lapply(split(pc_df, pc_df$pc), function(df){
    df$gene <- factor(df$gene, levels = df$gene)
    ggplot(df, aes(x = pc_loading, y = gene)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
      ggtitle(df$pc[1])
  })
  plot_grid_wrapper(gg_l, ncol = 8, align = 'v', axis = 'r', rel_height = 0.2
    , title = make_plot_title("Genes with highest PC loadings")
  )
  ggsave(paste0(out_graph, "pc_loading_plots.pdf"), width = 20, height = 7)

  # PC variance plot
  stdev <- seurat_obj@reductions$pca@stdev
  tibble(stdev = stdev, pc = c(1:length(stdev))) %>%
    ggplot(aes(x = pc, y = stdev)) +
      geom_point() +
      ggtitle(
        make_plot_title("PC variance")
        )
  ggsave(paste0(out_graph, "pc_variance_plot.pdf"), width = 5, height = 5)

}
################################################################################

### Correlation of PCs and metadata

plot_correlation_matrix_of_pcs_and_metadata <- function(seurat_obj){

  print("plot_correlation_matrix_of_pcs_and_metadata")

  # make tibble of metadata and pc scores
  metadata_to_plot <- c(
    "cell_ids", "age", "clinical_dx"
    , "library_id", "pmi_h", "rin_after_staining_and_lcm"
    , "sample_number", "sex", "weight_g", "number_umi", "number_genes")
    # , "mean_gc_content", "mean_cds_length")
  metadata_tb <-
    Embeddings(object = seurat_obj, reduction = "pca")[,1:10] %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell_ids") %>%
      as_tibble() %>%
      # clean names
      set_names(., { gsub("V", "PC", colnames(.)) }) %>%
      # add metadata
      inner_join(.
        , seurat_obj[[metadata_to_plot]]
        ) %>%
      select(-cell_ids) %>%
      mutate_if(is.character, as.factor)

  # calculate correlation of pcs and metadata
  cor_m <- matrix(nrow = ncol(metadata_tb), ncol = ncol(metadata_tb))
  colnames(cor_m) <- rownames(cor_m) <- colnames(metadata_tb)
  for(i in 1:nrow(cor_m)){
    for(j in 1:ncol(cor_m)){
      x <- metadata_tb %>% pull(rownames(cor_m)[i])
      y <- metadata_tb %>% pull(colnames(cor_m)[j])
      if (class(x) == "numeric" & class(y) == "numeric") {
        r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
      } else {
        lmout <- lm(y~x)
        r <- sqrt(abs(summary(lmout)$adj.r.squared))
      }
      cor_m[i,j] <- r
    }
  }
  # fill in NA values from dependent variable being factor
  for (j in 1:ncol(cor_m)) {
    if (all(is.na(cor_m[ ,j])) == TRUE) {
      cor_m[ ,j] <- cor_m[j, ]
    }
  }

  # plot
  # reverse column and order for plotting
  cor_m <- cor_m[ ,ncol(cor_m):1]
  ggcorrplot(cor_m, type = "full", lab = TRUE, show.diag = TRUE
    , title = make_plot_title("|Spearman's rho| correlation values")
  )
  ggsave(paste0(out_graph, "correlation_matrix.png"), width = 11, height = 11)

}
################################################################################

### tSNE and clustering using different PCs

plot_tsne_clustering_pcs_test <- function(){

  print("plot_tsne_clustering_pcs_test")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  plot_tsne_clustering_loom <- function(tsne_loom_path, seurat_cluster_col){
    tsne_gg <- filtered_loom[[tsne_loom_path]][,] %>% t() %>% as_tibble %>% rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cluster = filtered_loom[[seurat_cluster_col]][] %>% as.character) %>%
      ggplot(aes(x = tsne1, y = tsne2, color = cluster)) +
        geom_point(size = 0.05, alpha = 0.5) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        ggplot_set_theme_publication +
        ggtitle(paste0(tsne_loom_path, "\n", seurat_cluster_col))
    return(tsne_gg)
  }

  plot_tsne_clustering_loom("col_attrs/tsne_cell_embeddings_pc1to25", "col_attrs/cluster_pc1to25_res08")

  # organize plot_tsne_clustering_loom arguments with tibble
  # note cannot include () after plot_tsne_clustering_loom in pmap call
  tibble(
    tsne_loom_path = c(
      "col_attrs/tsne_cell_embeddings_pc1to25"
        , "col_attrs/tsne_cell_embeddings_pc1to50"
        , "col_attrs/tsne_cell_embeddings_pc1to75"
        , "col_attrs/tsne_cell_embeddings_pc1to100"
      )
    , seurat_cluster_col = c(
      "col_attrs/cluster_pc1to25_res08"
      , "col_attrs/cluster_pc1to50_res08"
      , "col_attrs/cluster_pc1to75_res08"
      , "col_attrs/cluster_pc1to100_res06"
    )) %>%
    # plot
    pmap(., plot_tsne_clustering_loom) %>%
    plot_grid_wrapper(plotlist = ., ncol = 2, rel_height = 0.2
      , align = 'v', axis = 'r'
      , title = paste0(script_name
        , "\n\nSeurat cluster and tSNE with different PCs used")
      )
    ggsave(paste0(out_graph, "pcs_test_tsne.png")
      , width = 12, height = 9)


  filtered_loom$close_all()
}
################################################################################

### Plot tSNE colored by clustering with different Seurat resolutions

plot_tsne_clustering_resolution_test <- function(seurat_obj){

  print("plot_tsne_clustering_resolution_test")

  cluster_col_names <- c("RNA_snn_res.0.4", "RNA_snn_res.0.5"
    , "RNA_snn_res.0.6" , "RNA_snn_res.0.7", "RNA_snn_res.0.8")

  gg_l <- map(cluster_col_names, function(variable){
    gg_tb <-
      Embeddings(object = seurat_obj, reduction = "tsne") %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(tsne1 = tSNE_1, tsne2 = tSNE_2) %>%
      inner_join(.
        , seurat_obj[[variable]] %>% rownames_to_column("cell_ids"))
    # gg_tb <- gg_l[[3]]
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$tsne1
      , dim_2 = gg_tb$tsne2
      , variable_value = gg_tb[[variable]]
      , title = variable
      , legend_title = variable
      , size = 0.5
      , alpha = 0.5
      , guide_size = 2
    )
  })
  # set list names as metadata to use for file names
  names(gg_l) <- sapply(gg_l, function(gg){gg$labels$title})
  dir.create(paste0(out_graph, "resolution_test_tsne"))
  lapply(names(gg_l), function(name){
    gg_l[[name]] +
      # geom_point(size = 0.75, alpha = 0.25) +
      ggtitle(make_plot_title(paste0(
        "Cells colored by ", name)))
    ggsave(paste0(out_graph, "resolution_test_tsne/", name, ".png")
      , width = 9, dpi = 400)
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(
      "Cells colored by cluster")
    )
  ggsave(paste0(out_graph, "resolution_test_tsne.png")
    , width = 28, height = 12, dpi = 200)

}
################################################################################

### dim red plot colored by metadata

plot_dim_reduction_colored_by_metadata <- function(
  seurat_obj, reduction = "tsne"){

  print("plot_dim_reduction_colored_by_metadata")

  metadata_to_plot <- c("cluster_pc1to100_res06"
    , "subclustering"
    , "sample_number", "age", "autopsy_id", "clinical_dx", "finalsite"
    , "library_id", "npdx1", "pmi_h", "prep", "rin_after_staining_and_lcm"
    , "region", "sex", "weight_g", "number_umi", "number_genes")
    # , "mean_gc_content", "mean_cds_length")
  metadata_to_plot <- metadata_to_plot[
    metadata_to_plot %in% colnames(seurat_obj[[]])]

  gg_l <- map(metadata_to_plot, function(variable){
    gg_tb <-
      Embeddings(object = seurat_obj, reduction = reduction) %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(dim_1 = 2, dim_2 = 3) %>%
      inner_join(.
        , seurat_obj[[variable]] %>% rownames_to_column("cell_ids"))
    # gg_tb <- gg_l[[3]]
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$dim_1
      , dim_2 = gg_tb$dim_2
      , variable_value = gg_tb[[variable]]
      , title = variable
      , legend_title = variable
      , size = 0.5
      , alpha = 0.5
      , guide_size = 2
    )
  })
  # set list names as metadata to use for file names
  names(gg_l) <- sapply(gg_l, function(gg){gg$labels$title})
  dir.create(paste0(out_graph, "metadata_", reduction))
  lapply(names(gg_l), function(name){
    gg_l[[name]] +
      # geom_point(size = 0.75, alpha = 0.25) +
      ggtitle(make_plot_title(paste0(
        "Cells colored by ", name
        , "\nDimensionality reduction: ", reduction)))
    ggsave(paste0(out_graph, "metadata_", reduction, "/", name, ".png")
      , width = 9, dpi = 400)
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0(
      "Cells colored by metadata"
      , "\nDimensionality reduction: ", reduction))
    )
  ggsave(paste0(out_graph, "metadata_", reduction, ".png")
    , width = 28, height = 22, dpi = 200)
    # , width = 6, height = 6, dpi = 200)
}
################################################################################

### dimensionality reduction colored by cluster

plot_dim_reduction_colored_by_cluster <- function(
  seurat_obj, cluster_col_name = "cluster_pc1to100_res06", reduction = "tsne"){

  print("plot_dim_reduction_colored_by_cluster")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  cluster_ids_tsne_tb <-
    Embeddings(object = seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    rownames_to_column("cell_ids") %>%
    as_tibble %>%
    rename(dim_1 = 2, dim_2 = 3) %>%
    inner_join(., seurat_obj[["cluster_ids"]] %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble() %>%
      mutate(cluster_ids = cluster_ids %>%
        as.character() %>%
        as.numeric())
    )
  max_cluster_number <- cluster_ids_tsne_tb$cluster_ids %>%
    max()
  gg_l <- map(0:max_cluster_number, function(cluster){
    cluster_ids_tsne_tb %>%
      mutate(cluster_membership = if_else(
        cluster_ids == cluster, TRUE, FALSE)) %>%
        # pull(cluster_membership) %>% table
    {plot_dim_reduction_colored_by_variable(
      dim_1 = .$dim_1
      , dim_2 = .$dim_2
      , variable_value = .$cluster_membership
      , title = paste0("Cluster: ", cluster)
      , legend_title = paste0("Cluster: ", cluster)
      , size = 0.1
      , alpha = 0.5
      , guide_size = 2
    ) + scale_color_manual(values = c("grey", "#00BFC4"))}
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = make_plot_title(paste0("Colored by cluster", "\n", reduction))
    )
  ggsave(paste0(out_graph, "cluster_", reduction, ".png")
    , width = 26, height = 4.25+max_cluster_number, dpi = 200)
    # , width = 6, height = 6, dpi = 200)
}
################################################################################

### Stacked bar charts of medadata by cluster

plot_cluster_metrics <- function(
  seurat_obj, cluster_col_name = "cluster_pc1to100_res06"){

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  metadata_tb <- seurat_obj[[]] %>%
    as_tibble %>%
    select_if(negate(is.numeric)) %>%
    # make sure clusters plot 0 to highest number
    mutate(cluster_ids = factor(
      cluster_ids
      , levels = cluster_ids %>%
        as.character() %>%
        as.numeric() %>%
        unique() %>%
        sort
      )) %>%
    gather(variable, value, -cell_ids, -cluster_ids, -orig.ident)

  metadata_tb %>%
    # mutate(cluster_ids = as.numeric(cluster_ids)) %>%
    filter(variable == "library_id") %>%
    group_by(cluster_ids) %>%
    summarise(count = n()) %>%
    ggplot(aes(x = cluster_ids, y = count)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = count, x = cluster_ids, y = count)
        , vjust = 0) +
      ggtitle(make_plot_title("Number of cells per cluster"))
      ggsave(paste0(out_graph, "number_cells_per_cluster_bargraph.pdf")
        , width = 11)

  # Count stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      ggplot(gg, aes(x = cluster_ids, fill = value)) +
        facet_wrap(~variable) +
        geom_bar() +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
     # +
        # guide = guide_legend(nrow=2)) +
      ggsave(paste0(out_graph, "metadata_by_cluster_stackedbar.png")
        , width = 24, height = 54, limitsize = FALSE)

  # Percent stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      gg <- gg %>%
        group_by(cluster_ids) %>%
        count(value) %>%
        mutate(percent = n/sum(n)) %>%
        mutate(variable = name) %>%
        ungroup() %>%
        mutate(cluster_ids = as.numeric(as.character(cluster_ids))) %>%
        mutate(cluster_ids = factor(cluster_ids, levels = sort(unique(cluster_ids))))
      ggplot(gg, aes(x = cluster_ids, y = percent, fill = value)) +
        facet_wrap(~variable) +
        geom_bar(stat = "identity") +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "r", ncol = 2, rel_height = 0.1
      , title = make_plot_title("Metadata by cluster"))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_stackedbar.png")
        , width = 24, height = 54, limitsize = FALSE)

  # heatmap of percent of cells from each library for each cluster
  name <- "library_id"
  metadata_tb %>%
    filter(variable %in% name) %>%
    group_by(cluster_ids) %>%
    count(value) %>%
    mutate(percent = n/sum(n)*100) %>%
    mutate(variable = name) %>%
    ungroup() %>%
    complete(cluster_ids, value, variable
      , fill = list(n = 0, percent = 0)) %>%
    ggplot(aes(x = cluster_ids, y = value, fill = percent)) +
      geom_tile() +
      geom_text(aes(label = paste0(round(percent, 1), "%  (", n, ")"))) +
      scale_fill_gradient(
        low = "white", high = "red", na.value = "grey", limits = c(0,25)) +
      ylab(name) +
      ggplot_set_theme_publication +
      ggtitle(make_plot_title(
        paste0("Percent of cells from each library for each cluster"
          , "\nEach column sums to 100%, () indicates number of cells")))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_heatmap.png")
        , width = 30, height = 7)

    # table of variance in percentage of cells for each cluster by disease by library
    metadata_tb %>%
      filter(variable %in% c("library_id", "clinical_dx")) %>%
      spread(key = variable, value = value) %>%
      group_by(cluster_ids) %>%
      count(library_id, clinical_dx) %>%
      mutate(percent = n/sum(n)*100) %>%
      group_by(cluster_ids, clinical_dx) %>%
      summarize(min(percent), max(percent), var(percent), sd(percent)) %>%
      write_csv(path = paste0(out_table, "percent_library_by_dx_cluster.csv"))

}
################################################################################

### tSNE colored by PC score

plot_tsne_colored_by_pc_score <- function(seurat_obj){

  print("plot_tsne_colored_by_pc_score")

  gg_tb <- map(1:15, function(pc){
    pc_variable_name <- paste0("pc_", pc)
    gg_tb <-
      # collect tsne values
      Embeddings(object = seurat_obj, reduction = "tsne") %>%
      as.data.frame() %>%
      rownames_to_column("cell_ids") %>%
      as_tibble %>%
      rename(tsne1 = tSNE_1, tsne2 = tSNE_2) %>%
      # collect pca scores
      inner_join(.
        , Embeddings(object = seurat_obj, reduction = "pca")[,pc] %>%
          enframe() %>%
          as_tibble() %>%
          rename(cell_ids = name, pc_scores = value)
        )
    plot_dim_reduction_colored_by_variable(
      dim_1 = gg_tb$tsne1
      , dim_2 = gg_tb$tsne2
      , variable_value = gg_tb$pc_scores
      , title = pc_variable_name
      , legend_title = pc_variable_name
    )
  }) %>%
  plot_grid_wrapper(plotlist = ., ncol = 4, rel_height = 0.1
    , title = make_plot_title(
      "Aggregated P1 and P2 samples colored by PC score")
    )
  ggsave(paste0(out_graph, "pc_score_tsne.png")
    , width = 14, height = 12, dpi = 200)

}
################################################################################

### tSNE colored by expression of known marker genes

plot_marker_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_pc1to100_res06"
  , reduction = "tsne"){

  print("plot_marker_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))

  # append source to cell type names to distinguish groupings
  groupings <- with(marker_genes_tb
    , paste0(marker_for, "_", clean_strings(source)))

  marker_genes_tb$groupings <- groupings

  # plot mean of groups
  plot_dim_reduction_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    , groupings = groupings
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , reduction = reduction
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 3
    , limit_low = 0
    , title = make_plot_title(
      "Mean normalized expression of cell type markers")
  )
  ggsave(paste0(
    out_graph, "markers_gene_expression_", reduction, ".png")
    , width = 18, height = 18, dpi = 200)

  # plot individually and output in sub directory
  dir.create(paste0(out_graph, "markers_individually_gene_expression_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(marker_genes_tb, marker_genes_tb$groupings)
  lapply(names(marker_genes_ltb), function(name){
    marker_genes_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
        plot_dim_reduction_colored_by_expression(
          genes = marker_genes_tb$gene_symbol
          , seurat_obj = seurat_obj
          , seurat_cluster_col = "cluster_ids"
          , expression_slot = "data"
          , reduction = reduction
          , expression_color_gradient = TRUE
          , ncol = 4
          , limit_high = 3
          , limit_low = 0
          , title = make_plot_title(paste0(name
            , "\nNormalized expression of cell type markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_individually_gene_expression_", reduction, "/", name, ".png")
      , width = 18, height = 2+0.75*length(marker_genes_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

}
################################################################################

### dimensionality reduction plot colored by expression of astrocyte markers

plot_marker_astrocyte_expression_dim_reduction <- function(
  seurat_obj, cluster_col_name = "cluster_pc1to100_res06", reduction = "tsne"){

  print("plot_marker_expression_dim_reduction")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # clean
  astro_markers_tb <- astro_markers_tb %>%
    filter(! is.na(gene_symbol)) %>%
    mutate(gene = gsub("\\*", "", .$gene_symbol)) %>% as.data.frame()

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = astro_markers_tb$gene
    , groupings = astro_markers_tb$marker_for
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , expression_color_gradient = TRUE
    , reduction = reduction
    , ncol = 4
    , limit_high = 2
    , limit_low = 0
    , title = make_plot_title("Mean normalized expression of astrocyte markers")
  )
  ggsave(paste0(out_graph, "markers_astrocytes_expression_", reduction, ".png")
    , width = 18, height = 8, dpi = 200)

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = astro_markers_tb$gene
    , groupings = astro_markers_tb$marker_for
    , seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , seurat_cluster_col = "cluster_ids"
    , reduction = reduction
    , legend_title = "expression\nz-score"
    , zscore = TRUE
    , ncol = 4
    , limit_high = 1
    , limit_low = -1
    , title = make_plot_title("Mean normalized expression of astrocyte markers")
  )
  ggsave(paste0(
    out_graph, "markers_astrocytes_expression_zscore_", reduction, ".png")
    , width = 18, height = 8, dpi = 200)

  # plot individually and output in sub directory
  dir.create(paste0(
    out_graph, "markers_astrocytes_individually_expression_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(astro_markers_tb, astro_markers_tb$marker_for)
  lapply(names(marker_genes_ltb), function(name){
    astro_markers_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
    plot_dim_reduction_colored_by_expression(
      genes = astro_markers_tb$gene
      , seurat_obj = seurat_obj
      , seurat_cluster_col = "cluster_ids"
      , expression_slot = "data"
      , reduction = reduction
      , expression_color_gradient = TRUE
      , ncol = 4
      , limit_high = 2
      , limit_low = 0
      , title = make_plot_title(paste0(name
        , "\nNormalized expression of astrocyte markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_astrocytes_individually_expression_", reduction, "/"
      , name, ".png")
      , width = 18, height = 2+0.75*length(astro_markers_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

  # plot individually and output in sub directory
  dir.create(paste0(
    out_graph, "markers_astrocytes_individually_expression_zscore_", reduction))
  # set list names as metadata to use for file names
  marker_genes_ltb <- split(astro_markers_tb, astro_markers_tb$marker_for)
  lapply(names(marker_genes_ltb), function(name){
    astro_markers_tb <- marker_genes_ltb[[name]]
    tryCatch(
      {
    plot_dim_reduction_colored_by_expression(
      genes = astro_markers_tb$gene
      , seurat_obj = seurat_obj
      , expression_slot = "scale.data"
      , seurat_cluster_col = "cluster_ids"
      , legend_title = "expression\nz-score"
      , reduction = reduction
      , zscore = TRUE
      , ncol = 4
      , limit_high = 1
      , limit_low = -1
      , title = make_plot_title(paste0(name
        , "\nNormalized expression of astrocyte markers")))
      },
        error = function(cond){
          return(NA)
      })
    ggsave(paste0(
      out_graph, "markers_astrocytes_individually_expression_zscore_", reduction, "/"
      , name, ".png")
      , width = 18, height = 2+0.75*length(astro_markers_tb$gene_symbol)
      , dpi = 200, limitsize = FALSE)
  })

}
################################################################################

### tSNE colored by expression of genes of interest

plot_genes_of_interest_expression_tsnes <- function(
  seurat_obj, cluster_col_name = "cluster_pc1to100_res06"){

  print("plot_genes_of_interest_expression_tsnes")

  genes <- c("NFKB1", "STAT3", "IFNAR1", "C3", "C4A", "C4B")

  seurat_obj$cluster_ids <- seurat_obj[[cluster_col_name]]

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = genes
    , seurat_obj = seurat_obj
    , seurat_cluster_col = "cluster_ids"
    , expression_slot = "data"
    , expression_color_gradient = TRUE
    , ncol = 4
    , limit_high = 2
    , limit_low = 0
    , title = make_plot_title("Mean normalized expression of genes of interest")
  )
  ggsave(paste0(out_graph, "genes_of_interest_expression_tsne.png")
    , width = 18, height = 6.5, dpi = 200)

  # plot
  plot_dim_reduction_colored_by_expression(
    genes = genes
    , seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , seurat_cluster_col = "cluster_ids"
    , legend_title = "expression\nz-score"
    , zscore = TRUE
    , ncol = 4
    , limit_high = 1
    , limit_low = -1
    , title = make_plot_title("Mean normalized expression of genes of interest")
  )
  ggsave(paste0(out_graph, "genes_of_interest_expression_zscore_tsne.png")
    , width = 18, height = 6.5, dpi = 200)

}
################################################################################

### top cluster enriched genes expression heatmap

plot_cluster_enriched_genes_heatmap <- function(
  seurat_obj, cluster_enriched_de_tb){

  print("plot_cluster_enriched_genes_heatmap")

  # gene_group_tb <- cluster_enriched_de_tb %>%
  #   filter(avg_logFC > 0) %>%
  #   group_by(cluster) %>%
  #   top_n(n = 20) %>%
  #   select(gene, cluster)
  gene_group_tb <- cluster_enriched_de_tb %>%
    filter(log2_fold_change > 0) %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = log2_fold_change) %>%
    select(gene, cluster) %>%
    rename(group = cluster)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$group
    , cluster_col = "cluster_pc1to100_res06"
    , plot_title = make_plot_title("Top enriched genes for each cluster")
    , limit_high = 4
  )
  ggsave(paste0(out_graph, "cluster_top_enriched_heatmap.png")
    , height = 30)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$group
    , cluster_col = "cluster_pc1to100_res06"
    , plot_title = make_plot_title("Top enriched genes for each cluster")
    , limit_low = -1
    , limit_high = 1
    , zscore = TRUE
  )
  ggsave(paste0(out_graph, "cluster_top_enriched_heatmap_zscore.png")
    , height = 30)

}
################################################################################

### top cluster expressed genes expression heatmap

plot_cluster_top_expressed_genes_heatmap <- function(
  seurat_obj, top_expressed_genes_tb){

  print("plot_cluster_enriched_genes_heatmap")

  gene_group_tb <- top_expressed_genes_tb %>%
    group_by(cluster) %>%
    top_n(n = 40, wt = mean_expression)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$cluster
    , cluster_col = "cluster_pc1to100_res06"
    , plot_title = make_plot_title("Top expressed genes for each cluster")
    , limit_high = 4
  )
  ggsave(paste0(out_graph, "cluster_top_expressed_heatmap.png")
    , height = 6*length(unique(top_expressed_genes_tb$cluster))
    , limitsize = FALSE)

  plot_mean_expression_genes_heatmap(
    seurat_obj = seurat_obj
    , expression_slot = "scale.data"
    , genes_to_plot = gene_group_tb$gene
    , genes_groupings = gene_group_tb$cluster
    , cluster_col = "cluster_pc1to100_res06"
    , plot_title = make_plot_title("Top expressed genes for each cluster")
    , limit_low = -1
    , limit_high = 1
    , zscore = TRUE
  )
  ggsave(paste0(out_graph, "cluster_top_expressed_heatmap_zscore.png")
    , height = 6*length(unique(top_expressed_genes_tb$cluster))
    , limitsize = FALSE)

}
################################################################################

### cluster enriched genes violin plots

plot_cluster_enriched_genes_violin_plots <- function(
  seurat_obj, cluster_enriched_de_tb){

  print("plot_cluster_enriched_genes_violin_plots")

  gene_group_tb <- cluster_enriched_de_tb %>%
    filter(log2_fold_change > 0) %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = log2_fold_change) %>%
    select(gene, cluster) %>%
    rename(group = cluster)

  # cell id cluster id key
  cellid_clusterid_tb <- Idents(seurat_obj) %>%
    enframe(name = "cell_id", value = "cluster")

  # expression z-scores
  mean_expr_zscores_tb <-
    # get expression data and subset to genes of interest
    GetAssayData(seurat_obj, slot = "data")[
      rownames(GetAssayData(seurat_obj, slot = "data")) %in%
        gene_group_tb$gene, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    gather("cell_id", "expression", -gene) %>%
    # add gene group info
    right_join(., gene_group_tb, by = c("gene" = "gene")) %>%
    # add sub-clusters
    left_join(., cellid_clusterid_tb)

  mean_expr_zscores_tb %>%
    mutate(gene = factor(gene, levels = gene_group_tb$gene)) %>%
    # group_by(group) %>%
    ggplot(aes(x = cluster, y = expression, fill = cluster)) +
      facet_wrap(group~gene, scales = "free") +
      geom_violin() +
      ggtitle(make_plot_title(
        "Top 5 cluster enriched genes for each cluster")) +
      ggplot_set_theme_publication
  ggsave(paste0(out_graph, "cluster_enriched_violin.png")
    , height = 2+2*length(unique(cluster_enriched_de_tb$cluster))
    , width = 10)

}
################################################################################

### Percentage of cell types

make_table_of_percent_of_cell_types <- function(){
  print("make_table_of_percent_of_cell_types")
  cluster_annot_tb <- bind_rows(
    # excitatory
    # did not include cluster 1, want to exclude sample P2_7 due to quality
    tibble(cluster = c(3,5,8,11,13,16,17), cell_type = "excitatory")
    , tibble(cluster = c(10,12,14), cell_type = "inhibitory")
    , tibble(cluster = c(2,18,22,24), cell_type = "astrocyte")
    , tibble(cluster = c(0,7,19,20), cell_type = "oligo")
    , tibble(cluster = c(9), cell_type = "micro")
    , tibble(cluster = c(15), cell_type = "endothelial")
    , tibble(cluster = c(16), cell_type = "opc")
  ) %>% mutate(cluster = as.character(cluster))

  filtered_loom$get.attribute.df(MARGIN = 2) %>%
    as_tibble %>%
    select(cluster = cluster_pc1to100_res06) %>%
    left_join(., cluster_annot_tb) %>%
    filter(cell_type != "NA") %>%
    group_by(cell_type) %>%
    summarise (n = n()) %>%
    mutate(percent = (n / sum(n)) * 100) %>%
    write.csv(file = paste0(out_table, "cell_type_percent.csv"))
}
################################################################################

### main function

main_function()
# main_function_subclustered_astrocytes()
################################################################################
