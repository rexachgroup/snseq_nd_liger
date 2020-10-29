# 2020-10-28 Lawrence Chen
# Import + clean 2020-08-16 pci_filtered_seurat_object.rdat.
# - Map sample preparation batch to sequencing batch.
# - Fix cellranger variables that were interpreted as strings.
# - save as rds without compression.
require(Seurat)
require(tidyverse)

## Command args to input to loop through args and run DE as array job
cmd_args <- commandArgs(trailingOnly = TRUE)
print(paste0("cmd_args: ", cmd_args))

## Inputs
in_seurat <-
  "../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"

## Outputs
out_table <- file.path(
  "../analysis/pci_import/")
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
meta <- nd_so@meta.data %>% as_tibble(rownames = "rownames")
meta_batch <- meta %>%
    mutate(library_batch = str_match(library_id, "(.+)_\\d+")[, 2],
           seq_batch = seq_batch_mapping[library_batch]) %>%
    glimpse

meta_numeric <- c(
    "Estimated.Number.of.Cells",
    "Mean.Reads.per.Cell",
    "Median.Genes.per.Cell",
    "Number.of.Reads",
    "Total.Genes.Detected",
    "Median.UMI.Counts.per.Cell"
)

meta_percentage <- c(
    "Valid.Barcodes",
    "Sequencing.Saturation",
    "Q30.Bases.in.Barcode",
    "Q30.Bases.in.RNA.Read",
    "Q30.Bases.in.Sample.Index",
    "Q30.Bases.in.UMI",
    "Reads.Mapped.to.Genome",
    "Reads.Mapped.Confidently.to.Genome",
    "Reads.Mapped.Confidently.to.Intergenic.Regions",
    "Reads.Mapped.Confidently.to.Intronic.Regions",
    "Reads.Mapped.Confidently.to.Exonic.Regions",
    "Reads.Mapped.Confidently.to.Transcriptome",
    "Reads.Mapped.Antisense.to.Gene",
    "Fraction.Reads.in.Cells"
)

parse_percent <- function(x) as.numeric(gsub("%", "", x, fixed = TRUE)) / 100

meta_parse_fix <- meta_batch %>%
    mutate(across(all_of(meta_numeric), parse_number)) %>%
    mutate(across(all_of(meta_percentage), parse_percent)) %>%
    glimpse

nd_so@meta.data <- meta_parse_fix %>%
    as.data.frame() %>%
    column_to_rownames("rownames")

saveRDS(nd_so, file.path(out_table, "pci_seurat.rds"), compress = FALSE)
