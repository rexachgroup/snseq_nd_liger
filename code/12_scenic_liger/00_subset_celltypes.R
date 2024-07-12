# Extract and subsample celltypes for SCENIC analysis.
# The initial code samples n = (15% of cell counts per library in each cell type).
# Couldn't replicate inital sampling of cells with the set seed, so use a predefined list of cells and genes to subset.

library(tidyverse)
library(Seurat)
library(SCopeLoomR)

set.seed(0xABCDEF)
SOBJ_FILE <- "../../analysis/seurat/20200816/pci_filtered_seurat_object.rdat"
OUT_DIR <- "../../analysis/12_liger_scenic/"
CELLS_FILE <- "../../resources/scenic/cells_presampled.csv"
GENES_FILE <- "../../resources/scenic/genes_presampled.csv"

target_regions <- c("calcarine")
target_celltype <- c("microglia")

main <- function() {
    dir.create(OUT_DIR)
    cells_tb <- read_csv(CELLS_FILE)
    genes_tb <- read_csv(GENES_FILE)
    load(SOBJ_FILE)
    sobj <- nd_so
    
    comb_tb <- expand_grid(regions = target_regions, cluster_cell_type = target_celltype)
    pwalk(comb_tb, function(...) {
        cr <- list(...)
        cells_s <- cells_tb %>% filter(region == cr$region, cluster_cell_type == cr$cluster_cell_type)
        genes_s <- genes_tb %>% filter(region == cr$region, cluster_cell_type == cr$cluster_cell_type)
        sobj_s <- subset(sobj, cells = cells_s$cells, features = genes_s$genes)
        out_expr <- str_glue("{OUT_DIR}/{cr$region}/{cr$cluster_cell_type}/expr.loom")
        out_genes <- str_glue("{OUT_DIR}/{cr$region}/{cr$cluster_cell_type}/expr_genes.tsv")
        dir.create(dirname(out_expr), recursive = T)
        build_loom(out_expr, GetAssayData(sobj_s, slot = "counts"))
        write.table(genes_s$genes, out_genes, row.names = F, col.names = F, quote = F)
    })
}

if (!interactive()) main()
