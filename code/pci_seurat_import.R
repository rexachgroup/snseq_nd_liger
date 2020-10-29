# 2020-10-28 Lawrence Chen
# Import + clean 2020-08-16 pci_filtered_seurat_object.rdat.
# - Map sample preparation batch to sequencing batch.
# - save as rds without compression.
require(Seurat)
require(tidyverse)

## Command args to input to loop through args and run DE as array job
cmd_args <- commandArgs(trailingOnly = TRUE)
print(paste0("cmd_args: ", cmd_args))

## Inputs
in_seurat <-
  "../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"

## Variables
date <- format(Sys.Date(), "%Y%m%d")

## Outputs
out_table <- file.path(
  "../analysis/pci_import/", date, "/tables/")
print(out_table)

# Make directories
dir.create(out_table, recursive = TRUE, showWarnings = FALSE)

seq_batch_mapping <- c(
    "P1" = "seq_P1",
    "P2" = "seq_P1",
    "P3" = "seq_P2",
    "P4" = "seq_P2",
    "P5" = "seq_P3"
)

load(in_seurat)
meta <- nd_so@meta.data %>% as_tibble()
meta_batch <- meta %>%
    mutate(library_batch = str_match(library_id, "(.+)_\\d+")[, 2],
           seq_batch = seq_batch_mapping[library_batch]) %>%
    glimpse

nd_so[["seq_batch"]] <- meta_batch$seq_batch
nd_so[["library_batch"]] <- meta_batch$library_batch

saveRDS(nd_so, file.path(out_table, "pci_seurat.rds"), compress = FALSE)
