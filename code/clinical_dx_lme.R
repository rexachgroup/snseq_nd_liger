# Damon Polioudakis
# 2020-05-05
# Run a linear mixed effects model comparing expression by dx for each cell type and region
# running with control as the reference
#   mean control = beta0 (control)
#   mean AD = beta0 (control) + beta AD
#   etc

# Must load modules:
#  module load R/3.6.0

# Sample qsub
#   qsub -N dx_lme -t 1-30 -tc 5 -l h_data=256G,h_rt=24:00:00,highp -m n qsub_r_script.sh -p clinical_dx_lme.R
###########################################################################

rm(list = ls())
set.seed(27)

require(Seurat)
require(tidyverse)
require(lme4)
require(lmerTest)
require(broom.mixed)
require(Matrix)
source("function_library.R")
source("ggplot_theme.R")
source("seurat_function_library.R")

sessionInfo()

## Command args to input to loop through args and run DE as array job
cmd_args <- commandArgs(trailingOnly = TRUE)
print(paste0("cmd_args: ", cmd_args))

## Inputs
in_seurat <-
  "../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"
marker_genes_tb <- read_csv(
  "../resources/cluster_markers_mixed_20181019.csv")
marker_genes_refined_tb <- read_csv(
  "../resources/20200703_cell_markers_refined.csv")
astro_markers_tb <- read_csv(
  "../resources/AstrocyteMarkers_Brie_20190110.csv")
in_exc_marker <- "../resources/excitatory_markers_20191023.csv"

## Variables
date <- format(Sys.Date(), "%Y%m%d")
script_name <- paste0("clinical_dx_lme.R ", date)
graph_subtitle <- "P1-5 C1-5 I1-5"

## Outputs
out_table <- paste0(
  "../analysis/clinical_dx_lme/", date, "/tables/clinical_dx_lme_")
print(out_table)

# Make directories
dir.create(dirname(out_table), recursive = TRUE)
###########################################################################

main_function <- function(){

  print("main_function")

  print("loading seurat object...")
  load(in_seurat)
  print("done loading seurat object...")

  # Idents(nd_so) <- nd_so[["RNA_snn_res.0.6"]]

  # nd_so <- assign_cell_type_to_cell(
  #   seurat_obj = nd_so, marker_genes_tb = marker_genes_refined_tb)
  # nd_so <- assign_cell_type_cluster(
  #   seurat_obj = nd_so, cluster_col_name = "RNA_snn_res.5")

  # make tibble of de function args to loop through with array job
  run_de_args_tb <- expand.grid(
      cell_type = c(
        "astrocyte",
        "endothelia",
        "ependymal",
        "excitatory",
        "inhibitory",
        "microglia",
        "oligodendrocyte",
        "opc",
        "pericyte",
        "t_cell"
      ),
      region = c(
        "preCG",
        "calcarine",
        "insula"
      )
    ) %>%
    as_tibble()
  print(paste0("Dimensions of run_de_args_tb: ",
    paste0(dim(run_de_args_tb), collapse = " ")))

  # assign args to variables using array job index
  cell_type <- run_de_args_tb$cell_type[[as.numeric(cmd_args[1])]]
  region <- run_de_args_tb$region[[as.numeric(cmd_args[1])]]

  # calculate DE with linear mixed effects model
  lme_out_tb <- run_lme_by_cell_type_and_region(
    seurat_obj = nd_so,
    region = region,
    cell_type = cell_type,
    # may need to add terms to design matrix in function
    model_design = "expression ~ clinical_dx + log_number_umi + (1 | library_id)",
    down_sample_cells = 10000)

  # filter the lme output and format into tibble
  lme_out_tb <- filter_and_format_lme_output(
    lme_out_tb = lme_out_tb,
    seurat_obj = nd_so,
    region = region,
    cell_type = cell_type,
    model_design = model_design,
    percent_detected_filter = -1
  )

  # save DE table to csv
  out_dx_de_csv <- paste0(
      out_table,
      region, "_",
      gsub("/", "_", cell_type), ".csv")
  print(paste0("out csv path: ", out_dx_de_csv))
  write_csv(lme_out_tb, path = out_dx_de_csv)

  print("end of main_function")
}
###########################################################################

run_lme_by_cell_type_and_region <- function(
  seurat_obj,
  cell_type = "excitatory",
  region = "calcarine",
  model_design = "expression ~ clinical_dx + log_number_umi + (1 | library_id)",
  down_sample_cells = NULL){

  print("run_lme_by_cell_type_and_region")

  print(cell_type)
  print(region)

  # subset seurat object to clinical_dx and cell type of interest
  # need to assign in global environment for seurat subset() function for some reason
  assign("cell_type_seurat_subset",
    cell_type,
    envir = .GlobalEnv)
  assign("region_seurat_subset",
    region,
    envir = .GlobalEnv)
  # cell_type and clinical_dx variables need to be global for some reason
  seurat_obj <- subset(x = seurat_obj,
    # subset = cluster_cell_type == cell_type)
    subset = cluster_cell_type == cell_type_seurat_subset)
  seurat_obj <- subset(x = seurat_obj,
    subset = region == region_seurat_subset)
  seurat_obj[[c("region", "clinical_dx", "cluster_cell_type")]] %>%
    table() %>%
    print()

  # subset to X number of cells to help with mem issues
  if(! is.null(down_sample_cells)){
    if(ncol(seurat_obj) > down_sample_cells){
      seurat_obj <- subset(
        seurat_obj, cells = sample(Cells(seurat_obj), down_sample_cells))
      print(paste0("after down-sampling to: ", down_sample_cells, " cells"))
      seurat_obj[[c("region", "clinical_dx", "cluster_cell_type")]] %>%
        table() %>%
        print()
    }
  }

  # Linear mixed effects model
  # expression as dependent variable
  expr_m <- GetAssayData(seurat_obj, slot = "data")
  # independent variables
  independent_vars_df <- data.frame(
    seurat_obj[[c("clinical_dx", "library_id", "number_umi"), drop = FALSE]])
  # log transform number_umi
  independent_vars_df$log_number_umi <- log(independent_vars_df$number_umi)
  # running with control as the reference
  #   mean control = beta0 (control)
  #   mean AD = beta0 (control) + beta AD
  #   etc
  independent_vars_df$clinical_dx <- independent_vars_df$clinical_dx %>%
    factor %>%
    relevel(., ref = "Control")
  # model design
  # remove library_chemistry for regions with only v2 or v3
  if(independent_vars_df$library_chemistry %>% unique() %>% length() < 2){
    mod <- gsub(" \\+ library_chemistry", "", model_design)
  } else {
    mod <- model_design
  }
  print(mod)
  # run
  lme_out_obj_l <- run_linear_mixed_effects_model(
    expr_m = expr_m, independent_vars_df = independent_vars_df, mod = mod)

  # format lme output into tibble
  lme_out_tb <-
    map_df(lme_out_obj_l, function(lme_out_obj){
      tidy(lme_out_obj, effects = "fixed")}, .id = "gene") %>%
    filter(grepl("dx", .$term) | is.na(term)) %>%
    select(-x, -effect) %>%
    filter(! is.na(term))

  print("end of run_lme_by_cell_type_and_region")

  return(lme_out_tb)
}
###########################################################################

filter_and_format_lme_output <- function(
  lme_out_tb = lme_out_tb,
  seurat_obj,
  clinical_dx = c("PSP-S", "Control"),
  cell_type = "excitatory",
  region = "calcarine",
  model_design = "expression ~ clinical_dx + log_number_umi + (1 | library_id)",
  percent_detected_filter = -1
  ){

  print("filter_and_format_lme_output")

  lme_out_tb <-
    lme_out_tb %>%
      pivot_wider(id_cols = gene, names_from = term,
          values_from = c(estimate, statistic, p.value))

  # format lm output into a tibble
  lme_out_tb <- tibble(
    gene = lme_out_tb %>% pull(gene),
    beta_ad = lme_out_tb %>% pull(estimate_clinical_dxAD),
    beta_ftd = lme_out_tb %>% pull(estimate_clinical_dxbvFTD),
    beta_psp = lme_out_tb %>% pull(`estimate_clinical_dxPSP-S`),
    t_statistic_ad = lme_out_tb %>% pull(statistic_clinical_dxAD),
    t_statistic_ftd = lme_out_tb %>% pull(statistic_clinical_dxbvFTD),
    t_statistic_psp = lme_out_tb %>% pull(`statistic_clinical_dxPSP-S`),
    pvalue_ad = lme_out_tb %>% pull(p.value_clinical_dxAD),
    pvalue_ftd = lme_out_tb %>% pull(p.value_clinical_dxbvFTD),
    pvalue_psp = lme_out_tb %>% pull(`p.value_clinical_dxPSP-S`)
  )

  # calculate percent of cells expressing gene
  # subset seurat object to region and cell type of interest
  # need to assign in global environment for seurat subset() function for some reason
  assign("cell_type_seurat_subset", cell_type, envir = .GlobalEnv)
  assign("region_seurat_subset", region, envir = .GlobalEnv)
  # cell_type and clinical_dx variables need to be global for some reason
  seurat_obj <- subset(x = seurat_obj,
    # subset = cluster_cell_type == cell_type_seurat_subset)
    subset = cluster_cell_type == cell_type_seurat_subset)
  seurat_obj <- subset(x = seurat_obj,
    subset = region == region_seurat_subset)
  expr_m <- GetAssayData(seurat_obj, slot = "data")
  # calculate percent of dx cells expressing gene
  lme_out_tb <- add_dx_percent_of_cells_expressing_gene(
    lme_out_tb,
    expr_m,
    seurat_obj,
    dx = "Control",
    percent_expressing_var_name = "percent_detected_ctrl")
  lme_out_tb <- add_dx_percent_of_cells_expressing_gene(
    lme_out_tb,
    expr_m,
    seurat_obj,
    dx = "AD",
    percent_expressing_var_name = "percent_detected_ad")
  lme_out_tb <- add_dx_percent_of_cells_expressing_gene(
    lme_out_tb,
    expr_m,
    seurat_obj,
    dx = "bvFTD",
    percent_expressing_var_name = "percent_detected_ftd")
  lme_out_tb <- add_dx_percent_of_cells_expressing_gene(
    lme_out_tb,
    expr_m,
    seurat_obj,
    dx = "PSP-S",
    percent_expressing_var_name = "percent_detected_psp")

  # filter to genes expressed in greater than X percent of cells
  print(paste0("number of genes before percent detected filter: ", nrow(lme_out_tb)))
  lme_out_tb <- lme_out_tb %>%
    filter(percent_detected_ad > percent_detected_filter | percent_detected_ftd > percent_detected_filter | percent_detected_psp > percent_detected_filter | percent_detected_ctrl > percent_detected_filter)
  print(paste0("number of genes after percent detected filter: ", nrow(lme_out_tb)))

  # add ensembl
  #   lme_out_tb$ensembl <- convert_gene_symbols_to_ensembl_ids(
  #     gene_symbols = lme_out_tb$gene)

  # fdr_pvalue correct
  lme_out_tb$fdr_pvalue_ad <- p.adjust(lme_out_tb$pvalue_ad, method = "BH")
  lme_out_tb$fdr_pvalue_ftd <- p.adjust(lme_out_tb$pvalue_ftd, method = "BH")
  lme_out_tb$fdr_pvalue_psp <- p.adjust(lme_out_tb$pvalue_psp, method = "BH")
  # check
  table(lme_out_tb$pvalue_ad < 0.05)
  table(lme_out_tb$fdr_pvalue_ad < 0.05)
  table(lme_out_tb$pvalue_ftd < 0.05)
  table(lme_out_tb$fdr_pvalue_ftd < 0.05)
  table(lme_out_tb$pvalue_psp < 0.05)
  table(lme_out_tb$fdr_pvalue_psp < 0.05)

  # formating
  lme_out_tb <- lme_out_tb %>%
    # order columns
    dplyr::select(
      gene,
      beta_ad,
      beta_ftd,
      beta_psp,
      t_statistic_ad,
      t_statistic_ftd,
      t_statistic_psp,
      fdr_pvalue_ad,
      fdr_pvalue_ftd,
      fdr_pvalue_psp,
      percent_detected_ctrl,
      percent_detected_ad,
      percent_detected_ftd,
      percent_detected_psp,
      pvalue_ad,
      pvalue_ftd,
      pvalue_psp) %>%
    # order by gene sym
    arrange(-desc(gene)) %>%
    # convert factors to characters
    mutate_if(is.factor, as.character) %>%
    as_tibble() %>%
    mutate(
      region = region,
      cell_type = cell_type)

  return(lme_out_tb)

  print("end of filter_and_format_lme_output")

}
###########################################################################

## Fxn: linear mixed effects model with lme4 and lmeTest packages
# independent_vars_df:
# ExpCondition RIN.y     Seq.PC1      Seq.PC2
#   CP         8.4       0.04792498   -0.090448567
#   CP         8.1       0.53502697   -0.287629654
#   CP         8.0       0.18824922   -0.155651102
#   VZ         8.4       0.02529722   -0.100858264
#   VZ         8.7       0.45139297   0.856908177
#   VZ         9.1       0.27861748   -0.248868277
# mod: "expression ~ clinical_dx + log_number_umi + (1 | library_id)"
run_linear_mixed_effects_model <- function(expr_m, independent_vars_df, mod) {

  print("run_linear_mixed_effects_model")

  # run lme
  lme_out_obj_l <- apply(expr_m, 1, function(expr) {
      mod <- as.formula(mod)
      dat <- data.frame(expression = expr, independent_vars_df)
      # use tryCatch to return NA when model can't be fit for a gene
      tryCatch({
        lmer(mod, data = dat)},
        error = function(e){NA})
    })

  print("end of run_linear_mixed_effects_model")

  return(lme_out_obj_l)

}

convert_gene_symbols_to_ensembl_ids <- function(
  gene_symbols,
  in_biomart_gene_info = "../../RNAseq_singlecellfetal/source/BiomaRt_Compile_GeneInfo_GRCh38_Ensembl87.csv"){
  # Converts vector of gene symbols to ensembl IDs
  # gene symbols must be character format
  # Need to load:
  biomart_df <- read.csv(in_biomart_gene_info, header = TRUE)
  print("convert_gene_symbols_to_ensembl_ids")
  idx <- match(gene_symbols, biomart_df$hgnc_symbol)
  ens <- biomart_df$ensembl_gene_id[idx]
  gene_symbols[! is.na(ens)] <- as.character(ens[! is.na(ens)])
  return(ens)
}

# calculate percent of control cells expressing gene
add_dx_percent_of_cells_expressing_gene <- function(
  lme_out_tb,
  expr_m,
  seurat_obj,
  dx = "Control",
  percent_expressing_var_name = "percent_detected_ctrl"){

  print("add_dx_percent_of_cells_expressing_gene")

  # subset to cells from dx
  idx <- seurat_obj[["clinical_dx"]][,1] == dx
  subset_expr_m <- expr_m[ ,idx]
  # calculate percent detected for each gene
  percent_detected <- (rowSums(subset_expr_m > 0) / ncol(subset_expr_m)) * 100
  is.na(percent_detected) <- 0
  percent_detected <- round(percent_detected, 1)
  percent_detected_df <- enframe(percent_detected,
    name = "gene", value = percent_expressing_var_name)
  # add to lme output table
  lme_out_tb <- left_join(lme_out_tb, percent_detected_df, by = "gene")

  print("end of add_dx_percent_of_cells_expressing_gene")

  return(lme_out_tb)
}
###########################################################################

main_function()
###########################################################################
