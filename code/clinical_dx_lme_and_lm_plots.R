# Damon Polioudakis
# 2020-07-18
# Plot down-sampled differential expression by clinical dx

# Must load modules:
#  module load R/3.6.0

# Sample qsub
#   qsub -N pl_dx_lm -t 1 -l h_data=128G,h_rt=12:00:00 -m n qsub_r_script.sh -p clinical_dx_lme_and_lm_plots.R
##########################################################################

set.seed(27)

require(Seurat)
require(cowplot)
require(viridis)
require(tidyverse)
require(readxl)
require(RColorBrewer)
require(ggpubr)
require(eulerr)
require(ComplexHeatmap)
require(patchwork)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")
sessionInfo()

## Inputs
in_lme_dir <- "../analysis/clinical_dx_lme/20200825/tables/"
in_lm_dir <- "../analysis/clinical_dx_lm/20200829/tables/"
#in_seurat <- "../analysis/seurat_pfc/PFC_mt8.rds"
in_seurat <- "../analysis/seurat/20200816/pci_filtered_seurat_object_test.rdat"
marker_genes_refined_tb <- read_csv(
  "../resources/20200703_cell_markers_refined.csv")
# disease genes curated by jessica
disease_genes_tb <- read_csv("../resources/disease_gene_markers_jessica_20200429.csv")
# genes curated by jessica
# selective_genes_tb <- read_csv("../resources/RNASEQ_SELECTIVE_GENES.csv")
selective_genes_tb <- read_xlsx("../resources/AGORAlist.xlsx")
# polo paper AD gwas genes
polo_ad_gwas_genes_tb <- read_csv("../resources/grubman_2019_st5_AD_GWAS_genes_ad.csv")
# mattis / tsai AD DE genes
in_tsai_de_dir <- "../resources/tsai_2019_de_genes/"

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- paste0("clinical_dx_lme_and_lm_plots.R ", date)
graph_subtitle <- "P1-5 C1-5 I1-5"
dx_order <- c("AD", "bvFTD", "PSP-S", "Control")
# dx_color <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
dx_color <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
region_order <- c("calcarine", "preCG", "insula")
region_color <- c("#ffffb3", "#fdb462", "#fccde5")

## Outputs
out_graph <- file.path("../",
  "/tmp/clinical_dx_lme_and_lm_")
out_table <- file.path("../",
  "/tmp/clinical_dx_lme_and_lm_")
print(out_graph)
print(out_table)

# Make directories
if (!dir.exists(dirname(out_graph))) {
  dir.create(dirname(out_graph), recursive = TRUE)
}
if (!dir.exists(dirname(out_table))) {
  dir.create(dirname(out_table), recursive = TRUE)
}
##########################################################################

main_function <- function(){

  print("main_function")

  print("loading seurat object...")
  load(in_seurat)
  #nd_so <- readRDS(in_seurat)
  print("done loading seurat object...")

  # load LME csvs
  lme_tb <- list.files(in_lme_dir, full.names = TRUE) %>%
    .[! grepl("keep_cell_ids", .)] %>%
    map(read_csv) %>%
    bind_rows()
  # load LM csvs
  lm_tb <- list.files(in_lm_dir, full.names = TRUE) %>%
    .[! grepl("keep_cell_ids", .)] %>%
    map(read_csv) %>%
    bind_rows()
  # combine
  lme_and_lm_tb <- inner_join(lme_tb, lm_tb, by = c("gene", "ensembl", "region", "cell_type"), suffix = c("_lme", "_lm"))
  # wrangle long
  lme_and_lm_long_tb <- wrangle_lme_and_lm_tb_to_long_format(lme_and_lm_tb)


  ## Plots

  # plot_dx_genes_expression_by_donor_heatmap_zscore
  lme_and_lm_tb %>%
    group_split(region, cell_type, .keep = TRUE) %>%
    map(function(data) {
      region <- unique(data$region)
      cell_type <- unique(data$cell_type)
      writeLines(str_glue("========== {region} {cell_type} ==============="))
      tryCatch(
        plot_dx_genes_expression_by_donor_heatmap_zscore(
          seurat_obj = nd_so,
          lme_and_lm_tb = data,
          cell_type_to_plot = cell_type,
          region_to_plot = region,
          show_row_names = FALSE,
          plot_width = 30,
          plot_height = 10
          ),
          error = function(e) {print(e)}
        )
    })
#   # calcarine
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#       plot_dx_genes_expression_by_donor_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "calcarine",
#         show_row_names = FALSE,
#         plot_width = 24,
#         plot_height = 9
#       )
#     }
#   # insula
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#     tryCatch(
#       plot_dx_genes_expression_by_donor_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "insula",
#         show_row_names = FALSE,
#         plot_width = 24,
#         plot_height = 9
#       ),
#       error = function(e) NA)
#     }
#   # preCG
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#     tryCatch(
#       plot_dx_genes_expression_by_donor_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "preCG",
#         show_row_names = FALSE,
#         plot_width = 24,
#         plot_height = 9
#       ),
#       error = function(e) NA)
#   }
# 
#   # plot_dx_genes_expression_heatmap_zscore
#   # calcarine
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#     tryCatch(
#       plot_dx_genes_expression_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "calcarine",
#         show_row_names = FALSE,
#         plot_width = 9,
#         plot_height = 9
#       ),
#       error = function(e) NA)
#   }
#   # insula
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#     tryCatch(
#       plot_dx_genes_expression_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "insula",
#         show_row_names = FALSE,
#         plot_width = 9,
#         plot_height = 9
#       ),
#       error = function(e) NA)
#   }
#   # preCG
#   for(cell_type_to_plot in unique(lme_and_lm_tb$cell_type)){
#     print(cell_type_to_plot)
#     tryCatch(
#       plot_dx_genes_expression_heatmap_zscore(
#         seurat_obj = nd_so,
#         lme_and_lm_tb = lme_and_lm_tb,
#         cell_type_to_plot = cell_type_to_plot,
#         region_to_plot = "preCG",
#         show_row_names = FALSE,
#         plot_width = 9,
#         plot_height = 9
#       ),
#       error = function(e) NA)
#     }

##### number of genes text on plot may be broken (add_cell_n_text)
  plot_number_de_genes_bar_plot(seurat_obj = nd_so, lme_and_lm_tb = lme_and_lm_tb, add_cell_n_text = TRUE)

##### broken
  # plot_de_dx_venn_diagrams(lme_tb = lme_tb)

  #   plt <- plot_genes_by_cell_type_bar_plot(
  #     genes = disease_genes_tb %>%
  #       filter(disease == "AD") %>%
  #       pull(gene),
  #     lme_and_lm_long_tb = lme_and_lm_long_tb)
  #   ggsave(paste0(out_graph, "de_genes_barplot_disease_ad_genes.pdf"),
  #     height = 50, width = 9, limitsize = FALSE)
  # 
  #   plt <- plot_genes_by_cell_type_bar_plot(
  #     genes = disease_genes_tb %>%
  #       filter(disease == "FTD") %>%
  #       pull(gene),
  #     lme_and_lm_long_tb = lme_and_lm_long_tb)
  #   ggsave(paste0(out_graph, "de_genes_barplot_disease_ftd_genes.pdf"),
  #     height = 140, width = 9, limitsize = FALSE)
  # 
  #   plt <- plot_genes_by_cell_type_bar_plot(
  #     genes = disease_genes_tb %>%
  #       filter(disease == "PSP") %>%
  #       pull(gene),
  #     lme_and_lm_long_tb = lme_and_lm_long_tb)
  #   ggsave(paste0(out_graph, "de_genes_barplot_disease_psp_genes.pdf"),
  #     height = 30, width = 9, limitsize = FALSE)

  selective_genes_tb <- selective_genes_tb %>%
    filter(hgnc_symbol %in% lme_and_lm_long_tb$gene) %>%
    arrange(hgnc_symbol)
  pdf(paste0(out_graph, "de_genes_barplot_selective_genes.pdf"), 
     height = 12, 
      width = length(unique(lme_and_lm_long_tb$region)) * 10)
    selective_genes_tb$hgnc_symbol %>%
      walk(function(gene) {
        writeLines(str_glue(" ======== {gene} ======== "))
        print(plot_genes_by_cell_type_bar_plot(
          genes = gene,
          lme_and_lm_long_tb = lme_and_lm_long_tb
        ))
      })
  dev.off()
  
  #percent_detected_ctrl
  #   plt <- plot_genes_by_cell_type_bar_plot(
  #     genes = selective_genes_tb %>%
  #       pull(1),
  #     lme_and_lm_long_tb = lme_and_lm_long_tb)
  #   ggsave(paste0(out_graph, "de_genes_barplot_selective_genes.pdf"), plt,
  #     height = nrow(selective_genes_tb) * 10, width = 30, limitsize = FALSE)
  # 
  plot_gene_expression_by_donor_jitter_box_plots(
    seurat_obj = nd_so,
    genes = disease_genes_tb %>%
      filter(disease == "FTD", gene %in% rownames(nd_so), !duplicated(gene)) %>%
      pull(gene),
    file_subname = "disease_ftd_genes_",
    graph_subtitle_2 = "genes from jessica"
  )

  plt <- plot_dim_reduction_colored_by_expression(
    seurat_obj = nd_so,
    genes = c("ACE2", "TMPRSS2"),
    seurat_cluster_col = "cluster_ids",
    reduction = "umap", 
    expression_slot = "data",
    expression_color_gradient = TRUE,
    zscore = FALSE,
    limit_high = 5
  )
  arePlotted <- names(which(apply(GetAssayData(nd_so, "counts")[c("ACE2", "TMPRSS2"), ], 2, sum) > 0))
  subAssay <- GetAssayData(nd_so, "counts")[,arePlotted]
  saveRDS(subAssay, paste0(out_graph, "subassay.rds"))
  ggsave(paste0(out_graph, "expression_tsne.pdf"), plot = plt, height = 12, width = 3 * 12)

  print("end of main_function")
}
##########################################################################

wrangle_lme_and_lm_tb_to_long_format <- function(lme_and_lm_tb){

  print("wrangle_lme_and_lm_tb_to_long_format()")

  lme_and_lm_long_tb <- full_join(
      # betas
      lme_and_lm_tb %>%
        select(ensembl, gene, region, cell_type, beta_ad_lme, beta_ftd_lme, beta_psp_lme) %>%
        gather(key = "clinical_dx", value = "beta_lme", -gene, -ensembl, -region, -cell_type) %>%
        # clean up dx
        mutate(clinical_dx = if_else(
          str_detect(clinical_dx, "ad"), "AD", if_else(
            str_detect(clinical_dx, "ftd"), "bvFTD", if_else(
              str_detect(clinical_dx, "psp"), "PSP-S", "NA")))),
      # t stats
      lme_and_lm_tb %>%
        select(ensembl, gene, region, cell_type, t_statistic_ad_lme, t_statistic_ftd_lme, t_statistic_psp_lme) %>%
        gather(key = "clinical_dx", value = "t_statistic_lme", -gene, -ensembl, -region, -cell_type) %>%
        # clean up dx
        mutate(clinical_dx = if_else(
          str_detect(clinical_dx, "ad"), "AD", if_else(
            str_detect(clinical_dx, "ftd"), "bvFTD", if_else(
              str_detect(clinical_dx, "psp"), "PSP-S", "NA")))),
      by = c("gene", "ensembl", "region", "cell_type", "clinical_dx")
    ) %>%
    full_join(.,
      # fdr p value
      lme_and_lm_tb %>%
        select(ensembl, gene, region, cell_type, fdr_pvalue_ad_lm, fdr_pvalue_ftd_lm, fdr_pvalue_psp_lm) %>%
        gather(key = "clinical_dx", value = "fdr_pvalue_lm", -gene, -ensembl, -region, -cell_type) %>%
        # clean up dx
        mutate(clinical_dx = if_else(
          str_detect(clinical_dx, "ad"), "AD", if_else(
            str_detect(clinical_dx, "ftd"), "bvFTD", if_else(
              str_detect(clinical_dx, "psp"), "PSP-S", "NA")))),
      by = c("gene", "ensembl", "region", "cell_type", "clinical_dx")
    )

  print("end of... wrangle_lme_and_lm_tb_to_long_format()")

  return(lme_and_lm_long_tb)
}
##########################################################################

plot_number_de_genes_bar_plot <- function(seurat_obj, lme_and_lm_tb, add_cell_n_text = FALSE){

  # count number of cells per clinical_dx, region, cell type
  cell_n_tb <- seurat_obj[[c("clinical_dx", "region", "cluster_cell_type")]] %>%
    group_by(clinical_dx, region, cluster_cell_type) %>%
    tally() %>%
    mutate(join_key = paste(clinical_dx, region, cluster_cell_type, sep = "_")) %>%
    ungroup() %>%
    select(., join_key, n) %>%
    rename(n_cells = n)

  # count number of DE genes that pass fold change and / or fdr pvalue filters
  number_signif_genes_tb <-
    lme_and_lm_tb %>%
      # ad
      mutate(ad_signif_up = if_else(
        abs(t_statistic_ad_lme) >= 2 &
        fdr_pvalue_ad_lm <= 0.05 &
        beta_ad_lme > 0.25, TRUE, FALSE)) %>%
      mutate(ad_signif_down = if_else(
        abs(t_statistic_ad_lme) >= 2 &
        fdr_pvalue_ad_lm <= 0.05 &
        beta_ad_lme < -0.25, TRUE, FALSE)) %>%
      # ftd
      mutate(ftd_signif_up = if_else(
        abs(t_statistic_ftd_lme) >= 2 &
        fdr_pvalue_ftd_lm <= 0.05 &
        beta_ftd_lme > 0.25, TRUE, FALSE)) %>%
      mutate(ftd_signif_down = if_else(
        abs(t_statistic_ftd_lme) >= 2 &
        fdr_pvalue_ftd_lm <= 0.05 &
        beta_ftd_lme < -0.25, TRUE, FALSE)) %>%
      # psp
      mutate(psp_signif_up = if_else(
        abs(t_statistic_psp_lme) >= 2 &
        fdr_pvalue_psp_lm <= 0.05 &
        beta_psp_lme > 0.25, TRUE, FALSE)) %>%
      mutate(psp_signif_down = if_else(
        abs(t_statistic_psp_lme) >= 2 &
        fdr_pvalue_psp_lm <= 0.05 &
        beta_psp_lme < -0.25, TRUE, FALSE)) %>%
      select(region, cell_type, ad_signif_up, ad_signif_down, ftd_signif_up, ftd_signif_down, psp_signif_up, psp_signif_down) %>%
      gather(key = "dx_signif", value = "signif", -region, -cell_type) %>%
      group_by(region, cell_type, dx_signif) %>%
      count(signif, name = "number_of_genes") %>%
      filter(signif == TRUE) %>%
      # add dx
      mutate(clinical_dx = if_else(
        str_detect(dx_signif, "ad"), "AD", if_else(
          str_detect(dx_signif, "ftd"), "bvFTD", if_else(
            str_detect(dx_signif, "psp"), "PSP-S", "NA"))))

  # add cell type numbers to number DE genes table
  number_signif_genes_tb <- number_signif_genes_tb %>%
    mutate(join_key = paste(
      gsub("_Control", "", clinical_dx),
      region,
      cell_type,
      sep = "_")) %>%
    left_join(., cell_n_tb, by = "join_key")

  # make down genes negative counts
  number_signif_genes_tb <- number_signif_genes_tb %>%
    mutate(number_of_genes = if_else(str_detect(dx_signif, "down"), -1 * number_of_genes, 1 * number_of_genes))

  # plot
  gg <- ggplot(number_signif_genes_tb, aes(x = cell_type, y = number_of_genes, group = dx_signif)) +
    facet_wrap(~region, ncol = 3) +
    geom_bar(stat = "identity", position = position_dodge(), aes(fill = clinical_dx)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    { if (add_cell_n_text == TRUE){
      geom_text(position = position_dodge(width = 1), aes(label = n_cells, hjust = 0), angle = 90)
    }} +
    # scale_x_continuous(breaks = dat$number_of_cells) +
    coord_cartesian(ylim = c(
      min(number_signif_genes_tb$number_of_genes)*1.1,
      max(number_signif_genes_tb$number_of_genes)*1.1)) +
    ggtitle(make_plot_title("Number of DE genes"))
  #ggsave(paste0(out_graph, "number_de_genes_bar_plot.pdf"), width = 12, height = 7)
  return(gg)

}
##########################################################################

venn_diagram <- function(list1, list2, list3, labels, title,
  fill_colors = c("#e41a1c", "#377eb8", "#4daf4a")) {
  genes <- unique(unlist(sapply(list(list1, list2, list3), as.character)))
  gene_m <- cbind(
    genes_in_list1 = genes %in% list1,
    genes_in_list2 = genes %in% list2,
    genes_in_list3 = genes %in% list3)
  colnames(gene_m) <- labels
  venn <- euler(gene_m)
  plot(venn, counts = TRUE, quantities = list(fontsize = 28),
    fill = fill_colors,
    labels = FALSE,
    # labels = list(fontsize = 28),
    main = list(label = title, fontsize = 12)
  )
}

plot_dx_up_intersection_venn_diagram <- function(cell_type, region){

  dx_1_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx > 0.2, fdr_pvalue < 0.5,
      clinical_dx == "AD_Control") %>%
      pull(gene)
  dx_2_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx > 0.2, fdr_pvalue < 0.5,
      clinical_dx == "bvFTD_Control") %>%
      pull(gene)
  dx_3_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx > 0.2, fdr_pvalue < 0.5,
      clinical_dx == "PSP-S_Control") %>%
      pull(gene)

  venn_diagram(
    list1 = dx_1_up_genes,
    list2 = dx_2_up_genes,
    list3 = dx_3_up_genes,
    labels = c("AD", "FTD", "PSP"),
    fill_colors = c("#F8766D", "#00BA38", "#619CFF"),
    title = paste0("DEG up in ", region, " ", cell_type)
  )

}

plot_dx_down_intersection_venn_diagram <- function(cell_type, region){

  dx_1_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx < -0.2, fdr_pvalue < 0.5,
      clinical_dx == "AD_Control") %>%
      pull(gene)
  dx_2_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx < -0.2, fdr_pvalue < 0.5,
      clinical_dx == "bvFTD_Control") %>%
      pull(gene)
  dx_3_up_genes <- lme_tb %>%
    filter(cell_type == cell_type, region == region,
      beta_clinical_dx < -0.2, fdr_pvalue < 0.5,
      clinical_dx == "PSP-S_Control") %>%
      pull(gene)

  venn_diagram(
    list1 = dx_1_up_genes,
    list2 = dx_2_up_genes,
    list3 = dx_3_up_genes,
    labels = c("AD", "FTD", "PSP"),
    fill_colors = c("#F8766D", "#00BA38", "#619CFF"),
    title = paste0("DEG down in ", region, " ", cell_type)
  )

}
##########################################################################

plot_genes_by_cell_type_bar_plot <- function(genes, lme_and_lm_long_tb){
  gg <- lme_and_lm_long_tb %>%
    filter(gene %in% genes) %>%
    mutate(fdr_pvalue_signif = symnum(
      fdr_pvalue_lm, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", " "))) %>%
    mutate(t_stat_filter = if_else(abs(t_statistic_lme) >= 2, ">=2", "")) %>%
    mutate(plot_annot_text = paste0(fdr_pvalue_signif, " ", t_stat_filter)) %>%
    ggplot(aes(y = beta_lme, x = cell_type, fill = clinical_dx)) +
      facet_wrap(gene~region, ncol = 3, scales = "free_x") +
      # geom_point(aes(size = percent_detected_dx, color = clinical_dx), position = position_dodge(0.3)) +
      geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.75) +
      geom_text(position = position_dodge(0.8),
        aes(label = plot_annot_text, vjust = "right", hjust = "center"),
        angle = 90, size = 5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(-0.5, 0.5)) +
      ggtitle(make_plot_title(paste0(
        "Clinical Dx vs Control DE by cell type and region",
        "\n**** <= 0.0001, *** <= 0.001, ** <= 0.01, * <= 0.05",
        "\n>=2 indicates that the abs(lme t statistic) >= 2")))
    return(gg)
}

plot_genes_by_cell_type_percentage <- function(genes, lme_and_lm_long_tb) {
  gg <- lme_and_lm_long_tb %>%
    filter(gene %in% genes) %>%
    mutate(fdr_pvalue_signif = symnum(
      fdr_pvalue_lm, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", " "))) %>%
    glimpse

}
##########################################################################

plot_gene_expression_by_donor_jitter_box_plots <- function(
  seurat_obj = nd_so,
  genes,
  file_subname = "",
  graph_subtitle_2 = "",
  plot_width = 12,
  plot_height = 10
  ){

  print("plot_gene_expression_by_donor_jitter_box_plots()")

  mean_expr_tb <-
    FetchData(seurat_obj, slot = "data",
      vars = c(genes, "library_id", "clinical_dx", "region", "cluster_cell_type")) %>%
      as_tibble() %>%
      filter(cluster_cell_type == "microglia", region == "preCG") %>%
      gather(key = "gene", value = "expression", -library_id, -clinical_dx, -region, -cluster_cell_type) %>%
      group_by(gene) %>%
      add_count(library_id, name = "cell_number") %>%
      group_by(library_id, gene, region, clinical_dx) %>%
      summarize(mean_expression = mean(expression), cell_number, .groups = "drop") %>%
      mutate(clinical_dx = factor(clinical_dx, levels = c("AD", "PSP-S", "bvFTD", "Control")))

  browser()

  my_t_test_comparisons <- list(c(1, 4), c(2, 4), c(3, 4))
  mean_dx_plot <- ggplot(mean_expr_tb, 
      aes(x = clinical_dx, y = mean_expression, color = clinical_dx, label = cell_number)) +
    facet_wrap(gene~region, scales = "free_y") +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = 0.5, position = position_jitter(seed = 1)) +
    geom_text(position = position_jitter(seed = 1)) +
    expand_limits(y = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(aes(group = clinical_dx)
      , comparisons = my_t_test_comparisons
      , method = "t.test"
      , p.adjust.method = "none"
      , label = "p.signif"
      # , label.y = max(dat$percent)
    ) +
    ggtitle(make_plot_title(paste0(
      "Gene expression in microglia",
      "\n", graph_subtitle_2)))
    # geom_jitter(position = position_dodge2(), alpha = 0.5)
    ggsave(paste0(
      out_graph, "mean_expression_by_sample_",
      file_subname, "jitter_boxplot.pdf"),
    mean_dx_plot, width = plot_width, height = plot_height)

#   mean_expr_zscore_tb <-
#     FetchData(seurat_obj, slot = "data",
#       vars = c(genes, "cell_ids", "library_id", "clinical_dx", "region", "cluster_cell_type")) %>%
#       as_tibble() %>%
#       filter(cluster_cell_type == "microglia") %>%
#       filter(region == "preCG") %>%
#       mutate_at(genes, scale) %>%
#       gather(key = "gene", value = "expression", -library_id, -clinical_dx, -region, -cluster_cell_type) %>%
#       # gather(key = "gene", value = "expression", -library_id, -clinical_dx, -region, -cluster_cell_type, -mean_mapt_expression, -mapt_ranking) %>%
#       add_count(library_id, name = "cell_number") %>%
#       group_by(library_id, gene, region, clinical_dx, cell_number) %>%
#       summarize(mean_expression = mean(expression)) %>%
#       ungroup() %>%
#       mutate(clinical_dx = factor(clinical_dx, levels = c("AD", "PSP-S", "bvFTD", "Control")))
# 
  mean_expr_tb %>%
    group_by(gene) %>%
    mutate(mean_expr_zscore = scale(mean_expression)) %>%
    ungroup() %>%
    ggplot(aes(x = clinical_dx, y = mean_expr_zscore, color = clinical_dx, label = cell_number)) +
      facet_wrap(gene~region, scales = "free_y") +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha = 0.5, position = position_jitter(seed = 1)) +
      geom_text(position = position_jitter(seed = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(make_plot_title(paste0(
        "Gene expression z-scores in microglia",
        "\n", graph_subtitle_2)))
      # geom_jitter(position = position_dodge2(), alpha = 0.5)
      ggsave(paste0(
        out_graph, "mean_expression_by_sample_zscore_",
        file_subname, "jitter_boxplot.pdf"),
      width = plot_width, height = plot_height)

  print("end of... plot_gene_expression_by_donor_jitter_box_plots()")

}
##########################################################################

plot_dx_genes_expression_by_donor_heatmap_zscore <- function(
  seurat_obj, lme_and_lm_tb, cell_type_to_plot, region_to_plot, show_row_names = FALSE,
  plot_width = 16, plot_height = 9
  ){

  print("plot_dx_genes_expression_by_donor_heatmap_zscore()")

  cell_type_enquo <- enquo(cell_type_to_plot)
  region_enquo <- enquo(region_to_plot)

  ##### may need to change to t_statistic_ad_lme_lme
  genes <-
    bind_rows(
      # ad
      filter(lme_and_lm_tb,
        abs(t_statistic_ad_lme) >= 2,
        fdr_pvalue_ad_lm <= 0.05,
        abs(beta_ad_lme) >= 0.25),
      # ftd
      filter(lme_and_lm_tb,
        abs(t_statistic_ftd_lme) >= 2,
        fdr_pvalue_ftd_lm <= 0.05,
        abs(beta_ftd_lme) >= 0.25),
      # psp
      filter(lme_and_lm_tb,
        abs(t_statistic_psp_lme) >= 2,
        fdr_pvalue_psp_lm <= 0.05,
        abs(beta_psp_lme) >= 0.25),
    ) %>%
    filter(region %in% !!region_enquo, cell_type %in% !!cell_type_enquo) %>%
    pull(gene) %>%
    unique()

  mean_expression_m <- FetchData(
    seurat_obj, vars = c("cluster_cell_type","clinical_dx", "region", "library_id", genes)) %>%
    filter(cluster_cell_type %in% !!cell_type_enquo) %>%
    select(-cluster_cell_type) %>%
    pivot_longer(cols = c(-clinical_dx, -library_id, -region)) %>%
    group_by(clinical_dx, library_id, region, name) %>%
    summarize(mean_expression = mean(value)) %>%
    pivot_wider(names_from = c("clinical_dx", "region", "library_id"), values_from = mean_expression, names_sep = ".") %>%
    column_to_rownames("name") %>%
    as.matrix() %>%
    t() %>%
    scale(center = TRUE, scale = TRUE) %>%
    t()

  # column annotation plots
  # annotation table
  column_annotation_df <-
    seurat_obj[[]] %>%
       count(library_id, region, clinical_dx, name = "number_of_cells") %>%
       mutate(key = paste(clinical_dx, region, library_id, sep = ".")) %>%
       # filter to expression matrix
       right_join(., tibble(key = colnames(mean_expression_m)), by = "key") %>%
       select(-library_id) %>%
       # set column order
       mutate(clinical_dx = factor(clinical_dx, levels = dx_order)) %>%
       arrange(region, clinical_dx) %>%
       column_to_rownames("key")
  # order columns of expression matrix to match
  idx <- match(rownames(column_annotation_df), colnames(mean_expression_m))
  mean_expression_m <- mean_expression_m[ ,idx]
  # make col annotation object
  dx_colors <- dx_color
  names(dx_colors) <- dx_order
  region_colors <- region_color
  names(region_colors) <- region_order
  column_annotation_obj <-
    HeatmapAnnotation(
     df = as.data.frame(select(column_annotation_df, -number_of_cells)),
     number_of_cells = anno_barplot(
       pull(column_annotation_df, number_of_cells)),
     col = list(
       percent_dx = circlize::colorRamp2(c(0, 100), c("white", "#006d2c")),
       clinical_dx = dx_colors,
       region = region_colors),
     border = FALSE)

  # row annotation plots
  # annotation table
  row_annotation_df <-
    lme_and_lm_tb %>%
      filter(region %in% !!region_enquo, cell_type %in% !!cell_type_enquo) %>%
  ##### may need to change to t_statistic_ad_lme_lme
      mutate(ad_signif = if_else(
        abs(t_statistic_ad_lme) >= 2 & fdr_pvalue_ad_lm <= 0.05 & abs(beta_ad_lme) >= 0.25, "YES", "NO")) %>%
      mutate(ftd_signif = if_else(
        abs(t_statistic_ftd_lme) >= 2 & fdr_pvalue_ftd_lm <= 0.05 & abs(beta_ftd_lme) >= 0.25, "YES", "NO")) %>%
      mutate(psp_signif = if_else(
        abs(t_statistic_psp_lme) >= 2 & fdr_pvalue_psp_lm <= 0.05 & abs(beta_psp_lme) >= 0.25, "YES", "NO")) %>%
      # filter(ad_signif == YES) %>% as.data.frame
      select(gene, ad_signif, ftd_signif, psp_signif) %>%
      # filter(! is.na(ad_signif) | ! is.na(ftd_signif) | ! is.na(psp_signif)) %>%
      # filter to expression matrix
      right_join(., tibble(gene = rownames(mean_expression_m)), by = "gene") %>%
      column_to_rownames("gene")
  row_annotation_obj <- rowAnnotation(
    df = row_annotation_df,
    col = list(
      ad_signif = c("YES" = dx_color[1], "NO" = "lightgray"),
      ftd_signif = c("YES" = dx_color[2], "NO" = "lightgray"),
      psp_signif = c("YES" = dx_color[3], "NO" = "lightgray")))

  # plot with ComplexHeatmap
  complex_heatmap_obj <- Heatmap(
    mean_expression_m,
    name = "normalized expression z-score",
    col = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"), space = "sRGB"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    row_title = "Genes passing FDR for a dx",
    column_title = "Samples",
    column_title_side = c("bottom"),
    top_annotation = column_annotation_obj,
    column_split = column_annotation_df[c("region", "clinical_dx")],
    right_annotation = row_annotation_obj,
    border = FALSE
    # row_order = row_order
    )
  # add title
  # convert heatmap to grob to plot with cowplot
  gb_heatmap <- grid.grabExpr(draw(complex_heatmap_obj))
  title <- make_plot_title(paste0(
    "Expression by sample of genes differentially expressed dx compared to control for ", region_to_plot, " ", cell_type_to_plot,
    "\nLME abs(t_statistic) >=2 & LME abs(beta) >= 0.25 & LM FDR p-value <= 0.05"))
  title <- ggdraw() + draw_label(title)
  rel_height <- 0.2
  plot_grid(title, gb_heatmap, ncol = 1, rel_heights = c(rel_height, 1))
  ggsave(paste0(out_graph, "expression_by_donor_heatmap_zscore_", region_to_plot, "_", cell_type_to_plot, ".png"),
   width = plot_width, height = plot_height)

  print("end of... plot_dx_genes_expression_by_donor_heatmap_zscore()")
}
##########################################################################

plot_dx_genes_expression_heatmap_zscore <- function(
  seurat_obj, lme_and_lm_tb, cell_type_to_plot, region_to_plot, show_row_names = FALSE,
  plot_width = 16, plot_height = 9
  ){

  print("plot_dx_genes_expression_heatmap_zscore()")
  print(out_graph)

  cell_type_enquo <- enquo(cell_type_to_plot)
  region_enquo <- enquo(region_to_plot)

  ##### may need to change to t_statistic_ad_lme_lme
  genes <-
    bind_rows(
      # ad
      filter(lme_and_lm_tb,
        abs(t_statistic_ad_lme) >= 2,
        fdr_pvalue_ad_lm <= 0.05,
        abs(beta_ad_lme) >= 0.25),
      # ftd
      filter(lme_and_lm_tb,
        abs(t_statistic_ftd_lme) >= 2,
        fdr_pvalue_ftd_lm <= 0.05,
        abs(beta_ftd_lme) >= 0.25),
      # psp
      filter(lme_and_lm_tb,
        abs(t_statistic_psp_lme) >= 2,
        fdr_pvalue_psp_lm <= 0.05,
        abs(beta_psp_lme) >= 0.25),
    ) %>%
    filter(region %in% !!region_enquo, cell_type %in% !!cell_type_enquo) %>%
    pull(gene) %>%
    unique()

  # calculate mean by donor, then mean of donor means to get the grouped mean
  mean_expression_m <- FetchData(
    seurat_obj, vars = c("cluster_cell_type","clinical_dx", "region", "library_id", genes)) %>%
    filter(cluster_cell_type %in% !!cell_type_enquo) %>%
    select(-cluster_cell_type) %>%
    pivot_longer(cols = c(-clinical_dx, -library_id, -region)) %>%
    # mean by donor
    group_by(clinical_dx, library_id, region, name) %>%
    summarize(mean_expression = mean(value), .groups = "drop") %>%
    # grouped mean
    group_by(clinical_dx, region, name) %>%
    summarize(mean_expression = mean(mean_expression), .groups = "drop") %>%
    pivot_wider(names_from = c("clinical_dx", "region"), values_from = mean_expression, names_sep = ".") %>%
    column_to_rownames("name") %>%
    as.matrix() %>%
    t() %>%
    scale(center = TRUE, scale = TRUE) %>%
    t()

  # column annotation plots
  # annotation table
  column_annotation_df <-
    seurat_obj[[]] %>%
       count(region, clinical_dx, name = "number_of_cells") %>%
       mutate(key = paste(clinical_dx, region, sep = ".")) %>%
       # filter to expression matrix
       right_join(., tibble(key = colnames(mean_expression_m)), by = "key") %>%
       # set column order
       mutate(clinical_dx = factor(clinical_dx, levels = dx_order)) %>%
       arrange(region, clinical_dx) %>%
       column_to_rownames("key")
   # order columns of expression matrix to match
   idx <- match(rownames(column_annotation_df), colnames(mean_expression_m))
   mean_expression_m <- mean_expression_m[ ,idx]
  # make col annotation object
  dx_colors <- dx_color
  names(dx_colors) <- dx_order
  region_colors <- region_color
  names(region_colors) <- region_order
  column_annotation_obj <-
    HeatmapAnnotation(
     df = as.data.frame(select(column_annotation_df, -number_of_cells)),
     number_of_cells = anno_barplot(
       pull(column_annotation_df, number_of_cells)),
     col = list(
       percent_dx = circlize::colorRamp2(c(0, 100), c("white", "#006d2c")),
       clinical_dx = dx_colors,
       region = region_colors),
     border = FALSE)

  # row annotation plots
  # annotation table
  row_annotation_df <-
    lme_and_lm_tb %>%
      filter(region %in% !!region_enquo, cell_type %in% !!cell_type_enquo) %>%
  ##### may need to change to t_statistic_ad_lme_lme
      mutate(ad_signif = if_else(
        abs(t_statistic_ad_lme) >= 2 & fdr_pvalue_ad_lm <= 0.05 & abs(beta_ad_lme) >= 0.25, "TRUE", "NA")) %>%
      mutate(ftd_signif = if_else(
        abs(t_statistic_ftd_lme) >= 2 & fdr_pvalue_ftd_lm <= 0.05 & abs(beta_ftd_lme) >= 0.25, "TRUE", "NA")) %>%
      mutate(psp_signif = if_else(
        abs(t_statistic_psp_lme) >= 2 & fdr_pvalue_psp_lm <= 0.05 & abs(beta_psp_lme) >= 0.25, "TRUE", "NA")) %>%
      # filter(ad_signif == TRUE) %>% as.data.frame
      select(gene, ad_signif, ftd_signif, psp_signif) %>%
      # filter(! is.na(ad_signif) | ! is.na(ftd_signif) | ! is.na(psp_signif)) %>%
      # filter to expression matrix
      right_join(., tibble(gene = rownames(mean_expression_m)), by = "gene") %>%
      column_to_rownames("gene")
  row_annotation_obj <- rowAnnotation(
    df = row_annotation_df,
    col = list(
      ad_signif = c("TRUE" = dx_color[1], "NA" = "lightgray"),
      ftd_signif = c("TRUE" = dx_color[2], "NA" = "lightgray"),
      psp_signif = c("TRUE" = dx_color[3], "NA" = "lightgray")))

  # plot with ComplexHeatmap
  complex_heatmap_obj <- Heatmap(
    mean_expression_m,
    name = "mean expression z-score",
    col = circlize::colorRamp2(c(-3, 0, 3), c("blue", "white", "red"), space = "sRGB"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = show_row_names,
    show_column_names = TRUE,
    row_title = "Genes passing FDR for a dx",
    column_title = "Dx and region",
    column_title_side = c("bottom"),
    top_annotation = column_annotation_obj,
    column_split = column_annotation_df[c("region", "clinical_dx")],
    right_annotation = row_annotation_obj,
    border = FALSE
    # row_order = row_order
    )
  # add title
  # convert heatmap to grob to plot with cowplot
  gb_heatmap <- grid.grabExpr(draw(complex_heatmap_obj))
  title <- make_plot_title(paste0(
    "Expression by dx of genes differentially expressed dx compared to control for ", region_to_plot, " ", cell_type_to_plot,
    "\nLME abs(t_statistic) >=2 & LME abs(beta) >= 0.25 & LM FDR p-value <= 0.05",
    "\nGrouped mean: calculated mean expression by sample, then mean of sample means"))
  title <- ggdraw() + draw_label(title)
  rel_height <- 0.2
  plot_grid(title, gb_heatmap, ncol = 1, rel_heights = c(rel_height, 1))
  ggsave(paste0(out_graph, "expression_by_dx_heatmap_zscore_", region_to_plot, "_", cell_type_to_plot, ".png"),
   width = plot_width, height = plot_height)

  print("end of... plot_dx_genes_expression_heatmap_zscore()")
}
##########################################################################

main_function()
##########################################################################
