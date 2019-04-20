# Damon Polioudakis
# 2018-09-10
# Clustering and initial analysis of P1 dataset

rm(list = ls())
set.seed(27)
sessionInfo()

require(cellrangerRkit)
require(xlsx)
require(cowplot)
require(viridis)
source("Function_Library.R")
source("GGplot_Theme.R")

options(stringsAsFactors = FALSE)


## Inputs

# cell ranger
cellranger_pipestance_path <- "../data/20180907/P1"
p1_crmo <- load_cellranger_matrix(cellranger_pipestance_path)
p1_crar <- load_cellranger_analysis_results(cellranger_pipestance_path)
p1_de_df <- read.csv(paste0(cellranger_pipestance_path
  , "/outs/analysis/diffexp/graphclust/differential_expression.csv"))
# metadata
p1_mt_df <- read.xlsx("../metadata/PGC_group1_071018.xlsx", 1)


## Outputs

# Paths
out_graph <- "../analysis/graphs/P1_Analysis/P1_Analysis_"
out_table <- "../analysis/tables/P1_Analysis/P1_Analysis_"
out_data <- "../analysis/data/P1_Analysis/P1_Analysis_"

# Make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
dir.create(dirname(out_data), recursive = TRUE)


## Other variables

script_name <- "P1_Analysis.R"
################################################################################

### Format

format_metadata <- function(){
  p1_mt_df$"Sample.." <- as.character(p1_mt_df$"Sample..")
  p1_mt_df$Prep <- as.character(p1_mt_df$Prep)
  return(p1_mt_df)
}
p1_mt_df <- format_metadata()

add_meta_data <- function(crmo, metadata
  , metadata_match_col = NULL, col_name = NULL){
  print("add_meta_data")
  # browser()
  if(is.null(dim(metadata))){
    pData(p1_crmo)[[col_name]] <- metadata
  } else {
    sample_id <- gsub(".*-", "", pData(p1_crmo)$barcode)
    idx <- match(sample_id, metadata[ ,names(metadata) == metadata_match_col])
    metadata <- metadata[idx, ]
    for(i in 1:length(names(metadata))){
      col_name <- names(metadata)[i]
      pData(p1_crmo)[[col_name]] <- metadata[ ,col_name]
    }
  }
  return(p1_crmo)
}
p1_crmo <- add_meta_data(
  crmo = p1_crmo
  , metadata = p1_mt_df
  , metadata_match_col = "Sample.."
)
p1_crmo <- add_meta_data(
  crmo = p1_crmo
  , metadata = as.character(p1_crar$clustering$graphclust$Cluster)
  , col_name = "cluster"
)
################################################################################

### Heatmap of cluster enriched genes

cluster_result <- p1_crar$clustering$graphclust
# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(p1_crmo, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster
# DE determined by sSeq method (Yu, Huber, & Vitek, 2013) which employs a
# negative binomial exact test.
# prioritized_genes <- prioritize_top_genes(
#   p1_crmo, cluster_result$Cluster, "sseq", min_mean = 0.5)

prioritized_genes <- lapply(names(p1_de_df)[grep("Log2", names(p1_de_df))], function(cluster){
  subset_p1_de_df <- p1_de_df[ ,c("Gene.Name", cluster)]
  subset_p1_de_df <- subset_p1_de_df[order(-subset_p1_de_df[ ,cluster]), ]
  return(subset_p1_de_df$Gene.Name[1:10])
})
prioritized_genes <- unlist(prioritized_genes)[1:10]

# create values and axis annotations for pheatmap
gbm_pheatmap(log_gene_bc_matrix(p1_crmo), prioritized_genes, cells_to_plot, n_genes=3, limits=c(-1,2))
ggsave(paste0(out_graph, "top_10_genes_heatmap.png"))
################################################################################

### tSNE colored by metadata

plot_tsne_colored_by_feature <- function(tsne_1, tsne_2, feature, title = NULL){
  print("Plot_tSNE_Colored_By_Feature")
  gg_DF <- data.frame(tsne_1 = tsne_1, tsne_2 = tsne_2, feature = feature)
  # Plot
  gg <- ggplot(gg_DF, aes(x = tsne_1, y = tsne_2, col = feature)) +
    geom_point(size = 0.1, alpha = 0.5) +
    guides(colour = guide_legend(override.aes = list(size = 7))) +
    ggplot_set_theme_publication +
    ggtitle(title)
  if(class(gg_DF$feature) == "numeric"){gg <- gg + scale_color_viridis()}
  return(gg)
}

plot_tsne_colored_by_metadata <- function(){
  print("plot_tsne_colored_by_metadata")
  condition <- names(pData(p1_crmo))
  condition <- condition[! condition %in% c(
    "barcode", "Note", "NA.", "Prep", "Region")]
  gg_l <- lapply(condition, function(data_type){
    tryCatch(
      {
        plot_tsne_colored_by_feature(
          tsne_1 = p1_crar$tsne$"TSNE.1"
          , tsne_2 = p1_crar$tsne$"TSNE.2"
          # , feature = pData(p1_crmo)[ ,1]
          , feature = pData(p1_crmo)[ ,names(pData(p1_crmo)) == data_type]
          , title = data_type
        )
      }, error = function(condition){
        message(condition)
        return(NA)}
    )
  })
  plot_grid_wrapper(gg_l, ncol = 4, align = 'v', axis = 'r'
    , title = paste0(
      script_name
      , "\n\nAggregated P1 samples colored by metadata"))
  ggsave(paste0(out_graph, "metadata_tsne.png")
    , width = 22, height = 14, dpi = 200)
}

plot_tsne_colored_by_metadata()
################################################################################
