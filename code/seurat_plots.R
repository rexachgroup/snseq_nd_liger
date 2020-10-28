# Damon Polioudakis
# 2018-10-18
# Plots of seurat analysis

# Must load modules:
#   module load R/3.6.0

# submit as job:
# qsub -N srt_plot -l h_data=128G,h_rt=24:00:00,highp qsub_r_script.sh -p seurat_plots.R
###########################################################################

rm(list = ls())
set.seed(27)

require(methods)
require(Seurat)
require(Matrix)
require(cowplot)
require(ggcorrplot)
require(viridis)
require(RColorBrewer)
require(ggpubr)
require(tidyverse)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
sessionInfo()

## inputs
in_seurat <-
  "/u/home/d/dpolioud/project-geschwind/nucseq_nd/analysis/seurat/20200816/pci_filtered_seurat_object.rdat"
marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
marker_genes_refined_tb <- read_csv(
  "../resources/20200703_cell_markers_refined.csv")
astro_markers_tb <- read_csv(
  "../resources/AstrocyteMarkers_Brie_20190110.csv")
in_exc_marker <- "../resources/excitatory_markers_20191023.csv"
polo_ad_gwas_genes_tb <- read_csv("../resources/grubman_2019_st5_AD_GWAS_genes_ad.csv")

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- paste0("seurat_plots.R ", date)
graph_subtitle <- "P1-5 C1-5, 8% MT filter, regress out UMI"

## outputs
out_graph <- paste0("../analysis/seurat/", date, "/graphs/seurat_")
out_table <- paste0("../analysis/seurat/", date, "/tables/seurat_")

# make directories
dir.create(dirname(out_graph), recursive = TRUE)
dir.create(dirname(out_table), recursive = TRUE)
###########################################################################

### functions

main_function <- function(){

  print("main_function")

  print("loading seurat object...")
  load(in_seurat)
  print("done loading seurat object...")

  # for when Seurat v2.3 was used for running CCA (not yet in v3.0)
  # nd_so <- UpdateSeuratObject(nd_so)

  # set cluster IDs and cell IDs
  nd_so$cluster_ids <- nd_so[["RNA_snn_res.1"]]
  nd_so$cell_ids <- nd_so[[]] %>% rownames

  # set data types and factor orders
  nd_so[["prep"]] <- nd_so[["prep"]] %>% pull() %>%
    factor(., levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"))
  nd_so[["clinical_dx"]] <- nd_so[["clinical_dx"]] %>% pull() %>%
    factor(., levels = c("AD", "bvFTD", "PSP-S", "Control"))

  # # assign cell type to cells and clusters
  # nd_so <- assign_cell_type_to_cell(
  #   seurat_obj = nd_so, marker_genes_tb = marker_genes_refined_tb)
  # nd_so <- assign_cell_type_cluster(
  #   seurat_obj = nd_so, cluster_col_name = "cluster_ids")

  plot_percent_of_cell_types_by_clinical_dx_and_region_boxplot(
    seurat_obj = nd_so)

  # # plot_marker_astrocyte_expression_dim_reduction(
  # #   seurat_obj = nd_so, reduction = "umap")

  plot_genes_of_interest_expression_dim_reduction(
    seurat_obj = nd_so, reduction = "umap",
    genes = c("NFKB1", "STAT3", "IFNAR1", "C3", "C4A", "C4B", "SPI1"
     , "APOE", "MAPT"))

  # T cell and macrophage markers
  # T cells CD247 CD69
  # Macrophage CD163
  # CD103 = ITGAE
  plot_genes_of_interest_expression_dim_reduction(
   seurat_obj = nd_so, reduction = "umap", out_suffix = "immune_cell_markers",
   genes = c("CD247", "CD69", "CD163", "SIP1", "CD22L", "CCR7", "CD44", "ITGAE"))

  # hypothesis is that T cells might be altered more in FTD > AD > PSP and perivascular macrophage unregulated in AD>FTD> PSP
  # CD103 - marker of interferon gamma secreting CD8 T cells that drive persistent microglia activation!!
  # {can be seen with CD69, SIP1, CD22L, CCR7, CD44}
  # We should def see these up in the bvFTD PreCG!!
  plot_number_of_cells_expressing_genes_box_jitter_plot(
    seurat_obj = nd_so, genes = c("CD247", "CD69", "CD163", "CD103", "SIP1", "CD22L", "CCR7", "CD44"))

  plot_ad_gwas_genes_expression_dim_reduction(
    seurat_obj = nd_so, reduction = "umap")
  # # plot_cluster_enriched_genes_heatmap(
  #   # seurat_obj = nd_so, cluster_enriched_de_tb = cluster_enriched_tb)

}
###########################################################################

### main function

main_function()
# main_function_subclustered_astrocytes()
###########################################################################
