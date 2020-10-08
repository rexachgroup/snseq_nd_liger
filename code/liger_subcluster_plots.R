# Damon Polioudakis
# lchen
# Plots of seurat analysis
###########################################################################

rm(list = ls())
set.seed(0)

require(methods)
require(Seurat)
require(Matrix)
require(cowplot)
require(ggcorrplot)
require(ggpubr)
require(viridis)
require(RColorBrewer)
require(tidyverse)
require(ComplexHeatmap)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
sessionInfo()
# require(xlsx)

## Command line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)

## inputs
#in_seurat <- "../analysis/seurat_lchen/liger_subcluster_subset.rds"
in_seurat_metadata <- readRDS("../analysis/seurat_lchen/liger_subcluster_subset.rds")

marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
marker_genes_refined_tb <- read_csv(
  "../resources/20200703_cell_markers_refined.csv")
astro_markers_tb <- read_csv(
  "../resources/AstrocyteMarkers_Brie_20190110.csv")
in_exc_marker <- "../resources/excitatory_markers_20191023.csv"
in_inhibitory_marker <- "../resources/inhibitory_markers_20200711.csv"
in_dx_markers <- "../resources/disease_gene_markers_jessica_20200429.csv"
polo_ad_gwas_genes_tb <- read_csv("../resources/grubman_2019_st5_AD_GWAS_genes_ad.csv")


# make directories
#dir.create(dirname(out_graph), recursive = TRUE)
#dir.create(dirname(out_table), recursive = TRUE)
###########################################################################

### functions

main_function <- function(seurat_rdat){

  print("main_function")
  load(seurat_rdat)
  writeLines(seurat_rdat)

  ## The plotting functions expect these variables to be defined as globals, which is an
  ## "interesting" interaction when running main_function in a loop
  #date <<- format(Sys.Date(), "%Y%m%d")
  out_file <- basename(seurat_rdat) %>% gsub("\\..*", "", .)
  graph_subtitle = basename(out_file)

  global_exports <- list(
    "script_name" = "liger_subcluster_plots.R",
    "out_file" = out_file,
    "graph_subtitle" = graph_subtitle,
    "dx_order" = c("AD", "bvFTD", "PSP-S", "Control"),
    "out_graph" = paste0("../analysis/seurat_lchen/graphs/", out_file, "_"),
    "out_table" = paste0("../analysis/seurat_lchen/tables/", out_file, "_")
  )
  list2env(global_exports, envir = .GlobalEnv)

  # set data types and factor orders
  liger_so[["prep"]] <- liger_so[["prep"]] %>% pull() %>%
    factor(., levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
  liger_so[["clinical_dx"]] <- liger_so[["clinical_dx"]] %>% pull() %>%
    fct_relevel(dx_order)

  liger_so$cluster_ids <- liger_so[["liger_clusters"]]
  liger_so$cell_ids <- liger_so[[]] %>% rownames

  plot_dim_reduction_colored_by_metadata(
    seurat_obj = liger_so, reduction = "umap",
    plot_height = 8, plot_width = 24,
    metadata_to_plot = c(
      "liger_clusters", "library_id", "prep", "library_chemistry",
      "clinical_dx", "region")
    )

  plot_metadata_by_cluster_percent_stacked_bar_plot(
    seurat_obj = liger_so, cluster_col_name = "liger_clusters",
    plot_height = 16, plot_width = 16,
    metadata_to_plot = c(
      "liger_clusters", "library_id", "prep", "library_chemistry",
      "clinical_dx", "region", "cell_type", "sex"))

  # plot dx markers from Jessica
  dx_marker_genes_tb <- read_csv(in_dx_markers) %>%
    select(gene_symbol = gene, marker_for = disease)
  dx_marker_genes <- dx_marker_genes_tb %>% pull(gene_symbol)
  # umaps
  plot_genes_of_interest_expression_dim_reduction(
    seurat_obj = liger_so, genes = dx_marker_genes,
    reduction = "umap", cluster_col_name = "liger_clusters",
    plot_width = 18, plot_height = 38, rel_height = 0.1,
    out_suffix = "dx_markers"
    )
  # heatmap
  plot_heatmap_of_marker_expression_with_dx_annotation(
    seurat_obj = liger_so, cluster_col_name = "liger_clusters",
    marker_genes_tb = dx_marker_genes_tb,
    z_score = FALSE,
    row_title = "Dx marker genes (Jessica)",
    plot_title = "Expression of disease marker genes currated by Jessica",
    plot_width = 12, plot_height = 12,
    out_graph,
    out_graph_suffix = "dx_marker_heatmap.png")

  ## Cell type specific plots

  if(grepl("excitatory", seurat_rdat) == TRUE){

    plot_excitatory_marker_expression_dim_reduction(
      seurat_obj = liger_so, cluster_col_name = "liger_clusters",
      reduction = "umap", in_exc_marker = in_exc_marker)

    plot_heatmap_of_excitatory_marker_expression(
      seurat_obj = liger_so, cluster_col_name = "liger_clusters",
      out_graph,
      in_exc_marker = in_exc_marker)

  }

  if(grepl("inhibitory", seurat_rdat) == TRUE){

    marker_genes_tb <- read_csv(in_inhibitory_marker) %>% select(gene_symbol, marker_for)
    marker_genes <- marker_genes_tb %>% pull(gene_symbol)

    plot_genes_of_interest_expression_dim_reduction(
      seurat_obj = liger_so, genes = marker_genes,
      reduction = "umap", cluster_col_name = "liger_clusters"
      )

    plot_heatmap_of_marker_expression_with_dx_annotation(
      seurat_obj = liger_so, cluster_col_name = "liger_clusters",
      marker_genes_tb = marker_genes_tb,
      z_score = TRUE,
      row_title = "Marker genes (Inma / Lake et al.)",
      plot_title = "Expression of Inma / Lake et al. inhibitory subtype markers",
      plot_width = 12, plot_height = 8,
      out_graph = out_graph,
      out_graph_suffix = "inhibitory_marker_heatmap_zscore.png")

  }

}
###########################################################################

plot_heatmap_of_excitatory_marker_expression <- function(
  seurat_obj, cluster_col_name, in_exc_marker, out_graph){

  print("plot_heatmap_of_excitatory_marker_expression()")

  seurat_obj$liger_clusters <- seurat_obj[[cluster_col_name]]

  marker_genes_tb <- read_csv(in_exc_marker)
  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))
  marker_genes_tb <- marker_genes_tb %>% filter(source == "hodge")

  expr_m <- FetchData(seurat_obj, vars = c(marker_genes_tb$gene_symbol, "liger_clusters")) %>%
    as_tibble(rownames = "cell_id") %>%
    gather(key = "gene", value = "expression", -cell_id, -liger_clusters) %>%
    group_by(liger_clusters, gene) %>%
    summarise(mean_expression = mean(expression)) %>%
    pivot_wider(
      names_from = liger_clusters,
      values_from = mean_expression) %>%
    right_join(.,
      select(marker_genes_tb, gene = gene_symbol, marker_for), by = "gene") %>%
    arrange(marker_for) %>%
    mutate(gene = paste(gene, marker_for)) %>%
    select(-marker_for) %>%
    # remove NA rows (complex heatmap can't cluster with these present)
    filter_at(vars(-gene), all_vars(! is.na(.))) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    t() %>%
    scale(., center = TRUE, scale = TRUE) %>%
    t()

  # calculate percent of dx cells for each cell type
  column_annotation_df <-
    seurat_obj[[]] %>%
       count(liger_clusters, clinical_dx) %>%
       mutate(clinical_dx = paste0(clinical_dx, "_percent_cells")) %>%
       group_by(liger_clusters) %>%
       mutate(percent = (n / sum(n)) * 100) %>%
       select(-n) %>%
       pivot_wider(id_cols = liger_clusters, names_from = clinical_dx, values_from = percent) %>%
       ungroup() %>%
       mutate(number_of_cells = c(seurat_obj[["liger_clusters"]] %>%
       group_by(liger_clusters) %>%
       count() %>%
       pull(n))) %>%
       select(liger_clusters, AD_percent_cells, bvFTD_percent_cells, `PSP-S_percent_cells`, Control_percent_cells, number_of_cells) %>%
       column_to_rownames("liger_clusters")
   column_annotation_obj <-
     HeatmapAnnotation(
      percent_dx = as.matrix(select(column_annotation_df, -number_of_cells)),
      number_of_cells = anno_barplot(
        pull(column_annotation_df, number_of_cells)),
      col = list(percent_dx = circlize::colorRamp2(
        c(0, 100), c("white", "#006d2c"))),
      border = TRUE)

  # annotation table for heatmap color bars
  row_annotation_df <-
    marker_genes_tb %>%
      mutate(key = paste(gene_symbol, marker_for)) %>%
      filter(key %in% rownames(expr_m)) %>%
      arrange(marker_for) %>%
      column_to_rownames("key") %>%
      select(marker_for) %>%
      as.data.frame()
  colors <- colorRampPalette(
    brewer.pal(9, "Blues"))(length(row_annotation_df$marker_for %>% unique()))
  names(colors) <- row_annotation_df %>% arrange(marker_for) %>% pull(marker_for) %>% unique()
  row_annotation_obj <- rowAnnotation(
    df = row_annotation_df,
    col = list(marker_for = colors)
  )

  # plot with ComplexHeatmap
  complex_heatmap_obj <- Heatmap(
    expr_m,
    name = "normalized expression z-score",
    col = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_title = "Marker genes (Lake et al.)",
    column_title = "Liger clusters",
    column_title_side = c("bottom"),
    top_annotation = column_annotation_obj,
    right_annotation = row_annotation_obj,
    border = TRUE
    # row_order = row_order
    )
  # add title
  # convert heatmap to grob to plot with cowplot
  gb_heatmap <- grid.grabExpr(draw(complex_heatmap_obj))
  title <- make_plot_title("Expression of Lake et al. excitatory layer markers")
  title <- ggdraw() + draw_label(title)
  rel_height <- 0.2
  plot_grid(title, gb_heatmap, ncol = 1, rel_heights = c(rel_height, 1))
  ggsave(paste0(out_graph, "exc_layer_marker_heatmap_zscore.png"),
    width = 12, height = 9)

  print("end of... plot_heatmap_of_excitatory_marker_expression()")
}
###########################################################################

plot_heatmap_of_marker_expression_with_dx_annotation <- function(
  seurat_obj, cluster_col_name, marker_genes_tb,
  z_score = FALSE,
  row_title = "genes", plot_title = "",
  plot_width = 12, plot_height = 12, out_graph, out_graph_suffix){

  # marker_genes_tb
  #    gene_symbol marker_for
  #   <chr>       <chr>
  # 1 LHX6        mge_derived
  # 2 ADARB2      cge_derived

  print("plot_heatmap_of_marker_expression_with_dx_annotation()")

  seurat_obj$cluster <- seurat_obj[[cluster_col_name]]

  # clean up marker tibble
  marker_genes_tb <- marker_genes_tb %>%
    filter(! is.na(gene_symbol)) %>%
    filter(gene_symbol %in% rownames(seurat_obj))
  # combine genes that are markers for multiple
  duplicated_genes <- marker_genes_tb$gene_symbol[duplicated(marker_genes_tb$gene_symbol)]
  for(duplicated_gene in duplicated_genes){
    combined_marker_for_tb <- marker_genes_tb %>% filter(gene_symbol %in% duplicated_gene) %>%
      mutate(combined_marker_for = paste0(marker_for, collapse = "/")) %>%
      slice(1)
    marker_genes_tb <- marker_genes_tb %>%
      mutate(marker_for = if_else(gene_symbol == combined_marker_for_tb$gene_symbol, combined_marker_for_tb$combined_marker_for, marker_for))
  }
  marker_genes_tb <- distinct(marker_genes_tb)

  # get expression data
  expr_m <- FetchData(seurat_obj, vars = c(marker_genes_tb$gene_symbol, "cluster")) %>%
    as_tibble(rownames = "cell_id") %>%
    gather(key = "gene", value = "expression", -cell_id, -cluster) %>%
    group_by(cluster, gene) %>%
    summarise(mean_expression = mean(expression)) %>%
    pivot_wider(
      names_from = cluster,
      values_from = mean_expression) %>%
    right_join(.,
      select(marker_genes_tb, gene = gene_symbol, marker_for), by = "gene") %>%
    arrange(marker_for) %>%
    mutate(gene = paste(gene, marker_for)) %>%
    select(-marker_for) %>%
    # remove NA rows (complex heatmap can't cluster with these present)
    filter_at(vars(-gene), all_vars(! is.na(.))) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  if(z_score == TRUE){
    expr_m <- expr_m %>%
      t() %>%
      scale(., center = TRUE, scale = TRUE) %>%
      t()
    color_palette <- circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    legend_title <- "normalized expression z-score"
  } else {
    legend_title <- "normalized expression"
    color_palette <- circlize::colorRamp2(c(0, 3), c("white", "red"))
  }

  # column annotation plots
  # calculate percent of dx cells for each cell type
  column_annotation_df <-
    seurat_obj[[]] %>%
       count(cluster, clinical_dx) %>%
       mutate(clinical_dx = paste0(clinical_dx, "_percent_cells")) %>%
       group_by(cluster) %>%
       mutate(percent = (n / sum(n)) * 100) %>%
       select(-n) %>%
       pivot_wider(id_cols = cluster, names_from = clinical_dx, values_from = percent) %>%
       ungroup() %>%
       mutate(number_of_cells = c(seurat_obj[["cluster"]] %>%
       group_by(cluster) %>%
       count() %>%
       pull(n))) %>%
       select(cluster, AD_percent_cells, bvFTD_percent_cells, `PSP-S_percent_cells`, Control_percent_cells, number_of_cells) %>%
       column_to_rownames("cluster")
  column_annotation_obj <-
    HeatmapAnnotation(
     percent_dx = as.matrix(select(column_annotation_df, -number_of_cells)),
     number_of_cells = anno_barplot(
       pull(column_annotation_df, number_of_cells)),
     col = list(percent_dx = circlize::colorRamp2(
       c(0, 100), c("white", "#006d2c"))),
     border = TRUE)

  # row annotation plots
  row_annotation_df <-
    marker_genes_tb %>%
      mutate(key = paste(gene_symbol, marker_for)) %>%
      filter(key %in% rownames(expr_m)) %>%
      arrange(marker_for) %>%
      column_to_rownames("key") %>%
      select(marker_for) %>%
      as.data.frame()
  row_annotation_obj <- rowAnnotation(
    df = row_annotation_df)

  # plot with ComplexHeatmap
  complex_heatmap_obj <- Heatmap(
    expr_m,
    name = legend_title,
    col = color_palette,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_title = row_title,
    column_title = cluster_col_name,
    column_title_side = c("bottom"),
    top_annotation = column_annotation_obj,
    right_annotation = row_annotation_obj,
    border = TRUE
    # row_order = row_order
    )
  # add title
  # convert heatmap to grob to plot with cowplot
  gb_heatmap <- grid.grabExpr(draw(complex_heatmap_obj))
  title <- make_plot_title(plot_title)
  title <- ggdraw() + draw_label(title)
  rel_height <- 0.2
  plot_grid(title, gb_heatmap, ncol = 1, rel_heights = c(rel_height, 1))
  ggsave(paste0(out_graph, out_graph_suffix),
    width = plot_width, height = plot_height)

  print("end of... plot_heatmap_of_marker_expression_with_dx_annotation()")
}

###########################################################################

### main function
map(unique(in_seurat_metadata$filepath), function(seurat_rdat_path) {
  main_function(seurat_rdat_path)
})
###########################################################################
