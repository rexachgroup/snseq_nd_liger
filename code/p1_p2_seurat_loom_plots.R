# Damon Polioudakis
# 2018-10-18
# Plots of seurat analysis

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(methods)
require(Seurat)
require(Matrix)
require(cowplot)
require(cellrangerRkit)
require(loomR)
require(ggcorrplot)
require(viridis)
require(tidyverse)
source("function_library.R")
source("ggplot_theme.R")
# require(xlsx)


## inputs
raw_loom_path <- "/u/flashscratch/d/dpolioud/p1_p2_seurat/20181206/p1_p2_raw.loom"
filtered_loom_path <- "../analysis/p1_p2_seurat/20181206/analyzed_data/p1_p2_filtered.loom"
marker_genes_tb <- read_csv(
  "../data/resources/cluster_markers_mixed_20181019.csv")

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- "p1_p2_seurat.R"

## outputs
out_graph <- paste0("../analysis/p1_p2_seurat/", date, "/graphs/p1_p2_seurat_")
out_table <- paste0("../analysis/p1_p2_seurat/tables/", date, "/p1_p2_seurat_")
out_data <- paste0(
  "../analysis/p1_p2_seurat/analyzed_data/", date, "/p1_p2_seurat_")


# make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_data), recursive = TRUE)
################################################################################

### functions

main_function <- function(){
  print("main_function")
  # plot_qc_metrics()
  # plot_pca()
  # plot_correlation_matrix_of_pcs_and_metadata()
  plot_tsne_clustering_pcs_test()
  plot_tsne_colored_by_metadata()
  plot_tsne_colored_by_pc_score()
  plot_cluster_metrics()
  plot_marker_expression_tsnes()
  make_table_of_percent_of_cell_types()
}

clean_strings <- function(string_vector){
  print("clean_strings")
  cleaned <- string_vector %>%
    gsub("* ", "_", .) %>%
    gsub("\\.", "_", .) %>%
    gsub("\\(", "_", .) %>%
    gsub("\\)", "_", .) %>%
    gsub("\\+", "and", .) %>%
    gsub("#", "_number", .) %>%
    gsub("_$", "", .) %>%
    gsub("__", "_", .) %>%
    tolower
  return(cleaned)
}

plot_tsne_colored_by_variable <- function(
  tsne_1, tsne_2, variable_value, facet_variable
  , title = NULL, guide_size = 4, legend_title = NULL
  , alpha = 0.5, size = 0.1
  , expression_color_gradient = FALSE){

  print("plot_tsne_colored_by_variable")

  gg_tb <- tibble(tsne_1 = tsne_1, tsne_2 = tsne_2
    , value = variable_value)
  if (! missing(facet_variable)){
    gg_tb <- mutate(gg_tb, facet_variable = facet_variable)}
  # Plot
  if (class(gg_tb$value) %in% c("character", "factor")){
    gg <- ggplot(gg_tb, aes(
        x = tsne_1, y = tsne_2, shape = value, col = value, group = value)) +
      geom_point(size = size, alpha = alpha, fill = NA) +
      scale_shape_manual(name = legend_title, values = rep(0:6,100)) +
      # option to facet
      { if (! missing(facet_variable)){
        facet_wrap(~facet_variable, scales = "free")
      }} +
      guides(color = guide_legend(title = legend_title
        , override.aes = list(size = guide_size, alpha = 1)))
  } else {
    gg <- ggplot(gg_tb, aes(x = tsne_1, y = tsne_2)) +
        geom_point(size = size, alpha = alpha, aes(col = value)) +
      scale_shape_manual(values = rep(0:6,100)) +
      # option to facet
      { if (! missing(facet_variable)){
        facet_wrap(~facet_variable, scales = "free")
      }} +
      # color scale options
      { if (expression_color_gradient == TRUE){
        scale_color_gradientn(name = legend_title
          , colours = c(
            "lightgrey", "#fee090", "#fdae61", "#f46d43", "#ca0020")
          )
      } else {
        scale_color_viridis(name = legend_title)
      }}
  }
    gg <- gg +
      ggplot_set_theme_publication +
      ggtitle(title)

  return(gg)
}

plot_tsne_colored_by_expression <- function(
  genes
  , groupings
  , loom_file
  , tsne_loom_path
  , expression_loom_path
  , cluster_loom_path
  , title = NULL
  , limit_high = 1.5
  , limit_low = -1.5
  , ncol = 4){

  print("plot_tsne_colored_by_expression")

  # function to collect tsne and expression values
  collect_tsne_and_expression_values <- function(genes){
    # collect tsne values
    gg_tb <- loom_file[[tsne_loom_path]][,] %>%
      t() %>% as_tibble %>%
      rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cell_names = loom_file[["col_attrs/cell_names"]][]) %>%
      # add mean expression values to tibble
      inner_join(.
        # collect gene expression values and calculate mean
        , loom_file[[expression_loom_path]][
          , loom_file$row.attrs$gene_names[] %in% genes] %>%

          as_tibble %>%
          mutate(expression = rowMeans(.)) %>%
          mutate(expression = replace(
            expression, expression > limit_high, limit_high)) %>%
          mutate(expression = replace(
            expression, expression < limit_low, limit_low)) %>%
          mutate(cell_names = loom_file[["col_attrs/cell_names"]][])
        )
    return(gg_tb)
  }

  # collect tsne and expression values
  # if groupings argument exists then calculate mean expression for each group
  if (! missing(groupings)){
    tsne_expression_ltb <- lapply(split(genes, groupings), function(genes){
      collect_tsne_and_expression_values(genes = genes)
      })
    # order by groupings argument, with any missing groupings removed
    list_order <- unique(groupings)[
      unique(groupings) %in% names(tsne_expression_ltb)]
    tsne_expression_ltb <- tsne_expression_ltb[list_order]
  # else collect expression for each gene individualluy
  } else {
    tsne_expression_ltb <- lapply(genes, function(gene){
    collect_tsne_and_expression_values(genes = gene)
    })
  }

  # plot
  expression_gg_l <- lapply(names(tsne_expression_ltb), function(grouping){
    tsne_expression_tb <- tsne_expression_ltb[[grouping]]
    plot_tsne_colored_by_variable(
      tsne_1 = tsne_expression_tb$tsne1
      , tsne_2 = tsne_expression_tb$tsne2
      , variable_value = tsne_expression_tb$expression
      , title = grouping
      , legend_title = "expression"
      , alpha = 0.25
      , expression_color_gradient = TRUE
    )
  })

  # tsne colored by cluster
  cluster_tb <- loom_file[[tsne_loom_path]][,] %>%
    t() %>% as_tibble %>%
    rename(tsne1 = V1, tsne2 = V2) %>%
    mutate(cell_names = loom_file[["col_attrs/cell_names"]][]) %>%
    mutate(cluster = loom_file[[cluster_loom_path]][])
  cluster_gg <- plot_tsne_colored_by_variable(
    tsne_1 = cluster_tb$tsne1
    , tsne_2 = cluster_tb$tsne2
    , variable_value = cluster_tb$cluster
    # , title = grouping
    , legend_title = "cluster"
    , alpha = 0.25
    , guide_size = 2
  )

  # combine graphs and plot
  # extract legend - make sure it exists (NA genes plot with no legend)
  plot_legend <- lapply(expression_gg_l, function(x) tryCatch(
    get_legend(x), error = function(e) NA)) %>%
      .[! is.na(.)] %>% .[[1]]
  # remove legends
  expression_gg_l <- lapply(expression_gg_l, function(gg){
    gg + theme(legend.position = "none")})
  # remove cluster legends if there are more than 10
  if (length(unique(cluster_tb$cluster)) > 10){
    cluster_gg <- cluster_gg + theme(legend.position = "none")}
  # plot
  plot_grid_wrapper(append(list(cluster_gg), expression_gg_l)
  , align = 'v', axis = 'r', ncol = ncol, title = title) %>%
    # add legend
    plot_grid(., plot_legend, rel_widths = c(1, .1))
}
################################################################################

### plot QC metrics (histograms, density plots)

plot_qc_metrics <- function(){

  print("plot_qc_metrics")

  p1_raw_loom <- connect(
    filename = raw_loom_path
    , mode = "r+")

  filtered_loom <- connect(
    filename = filtered_loom_path
    , mode = "r+")

  ## Plots

  dat_tb <- p1_raw_loom$get.attribute.df(MARGIN = 2) %>% as_tibble %>%
    select(c("cell_names", "number_genes", "number_umi", "percent_mito"
    , "number_umi", "mean_gc_content", "mean_cds_length", "library_id")) %>%
    # rename(. ,replace = c("sample_number" = "sample_number")) %>%
    mutate(data_type = "raw data") %>%
      bind_rows(.
        , filtered_loom$get.attribute.df(MARGIN = 2) %>% as_tibble %>%
          select(
            c("cell_names", "number_genes", "number_umi", "percent_mito", "mean_gc_content", "mean_cds_length", "library_id")
          ) %>%
          mutate(data_type = "filtered data")
        ) %>%
      mutate(data_type = factor(
        data_type, levels = c("raw data", "filtered data")))
      # gather(key = metric, value = value, -data_type, -cell_names)

  # violin boxplots
  means_tb <- dat_tb %>%
  gather(variable, value, -data_type, -cell_names, -library_id) %>%
    group_by(data_type, variable) %>%
    summarize(
      value_mean = round(mean(value, na.rm = TRUE), 1)
      , value_median = round(median(value, na.rm = TRUE),1)
    ) %>%
    mutate(label = paste0(
      "mean:\n", value_mean, "\nmedian:\n", value_median))
  dat_tb %>%
    gather(variable, value, -data_type, -cell_names, -library_id) %>%
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
      ggtitle(paste(script_name, "\n\nQC metrics"))
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots of QC metrics"
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots: number of genes detected per cell by sample"
      , "\nx axis limits 0-6,000"
      ))
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots: number of UMI per cell by sample"
      , "\nx axis limits 0-10,000"
      ))
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots: percent MT per cell by sample"
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots: mean CDS length per cell by sample"
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
    , title = (paste0(script_name
      , "\n\nHistograms and density plots: mean gc content per cell by sample"
      ))
  )
  ggsave(paste0(out_graph, "qc_mean_gc_by_sample_histogram_density.png")
    , width = 22, height = 6)

  p1_raw_loom$close_all()
  filtered_loom$close_all()

}
################################################################################

### PCA plots

plot_pca <- function(){

  print("plot_pca")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  # Plot genes with highest PC loadings
  pc_df <- filtered_loom[["row_attrs/pca_gene_loadings"]][1:8,] %>%
    t() %>% as_tibble %>%
    mutate(gene = filtered_loom[["row_attrs/gene_names"]][]) %>%
    gather(., key = "pc", value = "pc_loading", -gene) %>%
    mutate(pc = gsub("V", "", pc)) %>%
    group_by(pc) %>%
    filter(! is.na(pc_loading)) %>%
    arrange(pc_loading) %>%
    filter(row_number() %in% c(1:5,(n() - 4):n())) %>%
    ungroup() %>% as.data.frame()
  gg_l <- lapply(split(pc_df, pc_df$pc), function(df){
    df$gene <- factor(df$gene, levels = df$gene)
    ggplot(df, aes(x = pc_loading, y = gene)) +
      geom_point() +
      xlab("Loading") +
      ylab("Gene") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
      ggtitle(paste0("PC", df$pc[1]))
  })
  plot_grid_wrapper(gg_l, ncol = 4, align = 'v', axis = 'r', rel_height = 0.2
    , title = paste0(script_name
      , "\n\nGenes with highest PC loadings")
  )
  ggsave(paste0(out_graph, "pc_loading_plots.pdf"), width = 11, height = 7)

  # PC variance plot
  DimElbowPlot(filtered_loom, reduction.type = "pca", dims.plot = 100) +
    ggtitle(
      paste0(script_name
        , "\n\nPC variance")
      )
  ggsave(paste0(out_graph, "pc_variance_plot.pdf"), width = 5, height = 5)

  filtered_loom$close_all()
}
################################################################################

### Correlation of PCs and metadata

plot_correlation_matrix_of_pcs_and_metadata <- function(){

  print("plot_correlation_matrix_of_pcs_and_metadata")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  # make tibble of metadata and pc scores
  metadata_to_plot <- c(
    "cell_names", "age", "clinical_dx"
    , "library_id", "pmi_h", "rin_after_staining_and_lcm"
    , "sample_number", "sex", "weight_g", "number_umi", "number_genes"
    , "mean_gc_content", "mean_cds_length")
  metadata_tb <- filtered_loom[["col_attrs/pca_cell_embeddings"]][1:10, ] %>%
    t() %>% as_tibble %>%
    # clean names
    set_names(., { gsub("V", "PC", colnames(.)) }) %>%
    add_column(cell_names = filtered_loom[["col_attrs/cell_names"]][]) %>%
    # add metadata
    inner_join(.
      , filtered_loom$get.attribute.df(MARGIN = 2) %>%
        select(c(metadata_to_plot))
      ) %>%
    select(-cell_names) %>%
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
    , title = paste0(script_name, "\n\n|Spearman's rho| correlation values")
  )
  ggsave(paste0(out_graph, "correlation_matrix.png"), width = 11, height = 11)

  filtered_loom$close_all()
}
################################################################################

### tSNE and clustering using different PCs

plot_tsne_clustering_pcs_test <- function(){

  print("plot_tsne_clustering_pcs_test")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  plot_tsne_clustering_loom <- function(tsne_loom_path, cluster_loom_path){
    tsne_gg <- filtered_loom[[tsne_loom_path]][,] %>% t() %>% as_tibble %>% rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cluster = filtered_loom[[cluster_loom_path]][] %>% as.character) %>%
      ggplot(aes(x = tsne1, y = tsne2, color = cluster)) +
        geom_point(size = 0.05, alpha = 0.5) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        ggplot_set_theme_publication +
        ggtitle(paste0(tsne_loom_path, "\n", cluster_loom_path))
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
    , cluster_loom_path = c(
      "col_attrs/cluster_pc1to25_res08"
      , "col_attrs/cluster_pc1to50_res08"
      , "col_attrs/cluster_pc1to75_res08"
      , "col_attrs/cluster_pc1to100_res08"
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

### tSNE colored by metadata

plot_tsne_colored_by_metadata <- function(){

  print("plot_tsne_colored_by_metadata")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  metadata_to_plot <- c("cluster_pc1to100_res08"
    , "sample_number", "age", "autopsy_id", "clinical_dx", "finalsite"
    , "library_id", "npdx1", "pmi_h", "prep", "rin_after_staining_and_lcm"
    , "region", "sex", "weight_g", "number_umi", "number_genes"
    , "mean_gc_content", "mean_cds_length")

  gg_l <- map(metadata_to_plot, function(variable){
    gg_tb <- filtered_loom[["col_attrs/tsne_cell_embeddings_pc1to100"]][,] %>%
      t() %>% as_tibble %>%
      rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cell_names = filtered_loom[["col_attrs/cell_names"]][]) %>%
      # mutate(cluster = filtered_loom[["col_attrs/cluster_pc1to100_res08"]][] %>%
      #   as.character) %>%
      inner_join(., filtered_loom$get.attribute.df(MARGIN = 2) %>%
        select(c("cell_names", variable)))
    plot_tsne_colored_by_variable(
      tsne_1 = gg_tb$tsne1
      , tsne_2 = gg_tb$tsne2
      , variable_value = gg_tb[[variable]]
      , title = variable
      , legend_title = variable
      , size = 0.1
      , alpha = 0.1
      , guide_size = 2
    )
  })
  # set list names as metadata to use for file names
  names(gg_l) <- sapply(gg_l, function(gg){gg$labels$title})
  dir.create(paste0(out_graph, "metadata_tsne"))
  lapply(names(gg_l), function(name){
    gg_l[[name]] +
      geom_point(size = 0.75, alpha = 0.25) +
      ggtitle(paste0(script_name
        , "\n\nAggregated samples colored by ", name))
    ggsave(paste0(out_graph, "metadata_tsne/", name, ".png")
      , width = 9, dpi = 400)
  })
  plot_grid_wrapper(plotlist = gg_l, ncol = 4, rel_height = 0.1
    , align = 'v', axis = 'r'
    , title = paste0(script_name
      , "\n\nAggregated P1 samples colored by metadata")
    )
  ggsave(paste0(out_graph, "metadata_tsne.png")
    , width = 26, height = 22, dpi = 600)
    # , width = 6, height = 6, dpi = 200)

  filtered_loom$close_all()
}
################################################################################

### Stacked bar charts of medadata by cluster

plot_cluster_metrics <- function(){

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  metadata_tb <- filtered_loom$get.attribute.df(MARGIN = 2) %>% as_tibble %>%
    select_if(is.character) %>%
    mutate(cluster_pc1to100_res08 = factor(
      cluster_pc1to100_res08
      , levels = cluster_pc1to100_res08 %>% as.numeric %>% unique %>% sort)) %>%
    gather(variable, value, -cell_names, -cluster_pc1to100_res08)

  metadata_tb %>%
    # mutate(cluster_pc1to100_res08 = as.numeric(cluster_pc1to100_res08)) %>%
    filter(variable == "library_id") %>%
    ggplot(aes(x = cluster_pc1to100_res08)) +
      geom_bar() +
      ggtitle(paste0(script_name, "\n\nNumber of cells per cluster"))
      ggsave(paste0(out_graph, "number_cells_per_cluster_bargraph.pdf")
        , width = 11)

  # Count stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      ggplot(gg, aes(x = cluster_pc1to100_res08, fill = value)) +
        facet_wrap(~variable) +
        geom_bar() +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "l", ncol = 2, rel_height = 0.05
      , title = paste0(script_name, "\n\nMetadata by cluster"))
     # +
        # guide = guide_legend(nrow=2)) +
      ggsave(paste0(out_graph, "metadata_by_cluster_stackedbar.png")
        , width = 24, height = 30)

  # Percent stacked bar chart
  metadata_tb %>%
    with(. ,split(., variable)) %>%
    {lapply(names(.), function(name){
      gg <- .[[name]]
      gg <- gg %>%
        group_by(cluster_pc1to100_res08) %>%
        count(value) %>%
        mutate(percent = n/sum(n)) %>%
        mutate(variable = name)
        mutate(cluster_pc1to100_res08 = as.numeric(cluster_pc1to100_res08))
      ggplot(gg, aes(x = cluster_pc1to100_res08, y = percent, fill = value)) +
        facet_wrap(~variable) +
        geom_bar(stat = "identity") +
        { if (length(unique(gg$value)) > 11){
          scale_fill_manual(name = name, values = colorRampPalette(
            brewer.pal(n = 11, name = "Set3"))(length(unique(gg$value))))
        } else {
          scale_fill_brewer(name = name, palette = "Set3")
        }}
    })} %>%
    plot_grid_wrapper(., align = "v", axis = "l", ncol = 2, rel_height = 0.05
      , title = paste0(script_name, "\n\nMetadata by cluster"))
      ggsave(paste0(out_graph, "metadata_by_cluster_percent_stackedbar.png")
        , width = 24, height = 30)

  filtered_loom$close_all()
}
################################################################################

### tSNE colored by PC score

plot_tsne_colored_by_pc_score <- function(){

  print("plot_tsne_colored_by_pc_score")

  filtered_loom <- connect(filename = filtered_loom_path
    , mode = "r+")

  gg_tb <- map(1:15, function(pc){
    pc_variable_name <- paste0("pc_", pc)
    gg_tb <- filtered_loom[["col_attrs/tsne_cell_embeddings_pc1to100"]][,] %>%
      t() %>% as_tibble %>%
      rename(tsne1 = V1, tsne2 = V2) %>%
      mutate(cell_names = filtered_loom[["col_attrs/cell_names"]][]) %>%
      inner_join(.
        , filtered_loom[["col_attrs/pca_cell_embeddings"]][pc, ,drop = FALSE]
          %>% t() %>% as_tibble %>%
          mutate(cell_names = filtered_loom[["col_attrs/cell_names"]][]) %>%
          rename(c("V1" = pc_variable_name))
        )
    plot_tsne_colored_by_variable(
      tsne_1 = gg_tb$tsne1
      , tsne_2 = gg_tb$tsne2
      , variable_value = gg_tb[[pc_variable_name]]
      , title = pc_variable_name
      , legend_title = pc_variable_name
    )
  }) %>%
  plot_grid_wrapper(plotlist = ., ncol = 4, rel_height = 0.1
    , title = paste0(script_name
      , "\n\nAggregated P1 samples colored by PC score")
    )
  ggsave(paste0(out_graph, "pc_score_tsne.png")
    , width = 14, height = 12, dpi = 200)

  filtered_loom$close_all()
}
################################################################################

### tSNE colored by expression of known marker genes

plot_marker_expression_tsnes <- function(){

  print("plot_marker_expression_tsnes")

  filtered_loom <- connect(filename = filtered_loom_path, mode = "r+")

  marker_genes_tb <- marker_genes_tb %>% filter(! is.na(gene_symbol))

  # append source to cell type names to distinguish groupings
  groupings <- with(marker_genes_tb
    , paste0(marker_for, "_", clean_strings(source)))

  # plot
  plot_tsne_colored_by_expression(
    genes = marker_genes_tb$gene_symbol
    , groupings = groupings
    , loom_file = filtered_loom
    , tsne_loom_path = "col_attrs/tsne_cell_embeddings_pc1to100"
    , expression_loom_path = "layers/norm_data"
    , cluster_loom_path = "col_attrs/cluster_pc1to100_res08"
    , ncol = 4
    , limit_high = 4
    , limit_low = 0
    , title = paste0(script_name
      , "\n\nMean normalized expression of cell type markers")
  )
  ggsave(paste0(out_graph, "marker_gene_expression_tsne.png")
    , width = 18, height = 18, dpi = 300)

  filtered_loom$close_all()
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
    select(cluster = cluster_pc1to100_res08) %>%
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
################################################################################
