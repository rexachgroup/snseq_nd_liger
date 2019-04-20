# Damon Polioudakis
# 2018-09-26
# Run Seurat on P1

# Must load modules:
#  module load gcc/4.9.3
#  module load R/3.3+
################################################################################

rm(list = ls())
set.seed(27)
sessionInfo()

require(Seurat)
require(loomR)
source("function_library.R")
source("ggplot_theme.R")

work_dir <- "/u/flashscratch/d/dpolioud/"

mca <- connect(filename = "mca.loom", mode = "r+")

download.file(
  "https://www.dropbox.com/s/8d8t4od38oojs6i/MCA.zip?dl=1"
  , destfile = paste0(work_dir, "mca.zip"))

dir.create(paste0(work_dir, "mca"), recursive = TRUE)

unzip(paste0(work_dir, "mca.zip")
  , exdir = paste0(work_dir, "mca")
  , files = NULL, list = FALSE, overwrite = TRUE, junkpaths = FALSE
  , unzip = "internal", setTimes = FALSE)


mca.matrix <- readRDS(file = paste0(work_dir, "mca/MCA/MCA_merged_mat.rds"))
mca.metadata <- read.csv(
  paste0(work_dir, "mca/MCA/MCA_All-batch-removed-assignments.csv")
  , row.names = 1)

# Only keep annotated cells
cells.use <- which(x = colnames(x = mca.matrix) %in% rownames(x = mca.metadata))
mca.matrix <- mca.matrix[, cells.use]
mca.metadata <- mca.metadata[colnames(x = mca.matrix), ]
# Create the loom file
mca <- create(filename = "mca.loom", data = mca.matrix
  , display.progress = FALSE, calc.numi = TRUE)
# Leaves us with 242k cells
mca

# Pull the tissue information for the cells we're analysing
tissues <- as.character(x = mca.metadata[, "Tissue"])
# Works similarly to AddMetaData for Seurat objects
mca$add.col.attribute(attribute = list(tissue = tissues))

methods(class = "loom")

NormalizeData(object = mca, chunk.size = 1000, scale.factor = 10000, display.progress = FALSE)

FindVariableGenes(object = mca)
hv.genes <- head(x = GetVariableGenes(object = mca)$index, n = 1000)

# Pull the indices of the mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = mca[["row_attrs/gene_names"]][], value = FALSE)
# Calculate percent.mito and store directly to the loom object Similar to
# creating a vector with the percentage of expression from mitochondrial
# genes, then using AddMetaData to put in to an object
mca$apply(name = "col_attrs/percent_mito", FUN = function(mat) {
    return(rowSums(x = mat[, mito.genes])/rowSums(x = mat))
}, MARGIN = 2, dataset.use = "matrix")
ScaleData(object = mca, genes.use = hv.genes, chunk.size = 20, display.progress = FALSE,
    vars.to.regress = "percent_mito")

RunPCA(object = mca, pc.genes = hv.genes, online.pca = FALSE, pcs.compute = 100,
    do.print = TRUE, pcs.print = 1:5, genes.print = 5, display.progress = FALSE)

PCElbowPlot(object = mca, num.pc = 100)
ggsave("../tmp/mca_pca_elbow_plot.png")

FindClusters(object = mca, reduction.type = "pca", dims.use = 1:75
  , resolution = 3, save.SNN = TRUE, n.start = 10, nn.eps = 0.5
  , print.output = FALSE)

RunTSNE(object = mca, reduction.use = "pca", dims.use = 1:75, max_iter = 2000, nthreads = 4, overwrite = TRUE)


p1 <- DimPlot(object = mca, reduction.use = "tsne", no.legend = TRUE, do.return = TRUE,
    pt.size = 0.1, ident.use = "col_attrs/res.3", vector.friendly = TRUE) +
    ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = mca, reduction.use = "tsne", no.legend = TRUE, do.return = TRUE,
    pt.size = 0.1, ident.use = "col_attrs/tissue", vector.friendly = TRUE) +
    ggtitle("Tissue") + theme(plot.title = element_text(hjust = 0.5))
plot_grid(p1, p2)
